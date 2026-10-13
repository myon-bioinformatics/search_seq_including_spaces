"""Bounded stdlib HTTP/webhook and MCP stdio adapters for deterministic tools.

No local file paths, URLs, network acquisition or model credentials are accepted
by tool calls. All protocols share the same explicit schemas and dispatcher.
"""
import argparse
from http.server import BaseHTTPRequestHandler, HTTPServer
import json
import math
import sys

from ambiguity_candidates import generate_candidates
from anchor_matcher import find_matches
from comparison import compare_aligned
from protein_ionization import ionization
from sequence_analysis import nucleotide_features, marker_clues, fold_rna
from structured_readers import strict_json

MAX_REQUEST = 1024 * 1024
MAX_RESPONSE = 2 * 1024 * 1024
VERSION = '2025-11-25'


def string(*, max_length=100000, enum=None):
    result={'type':'string','maxLength':max_length}
    if enum is not None: result['enum']=list(enum)
    return result


def integer(minimum=0,maximum=10000):
    return {'type':'integer','minimum':minimum,'maximum':maximum}


def definition(name, function, required, properties):
    return {'name':name,'description':function.__doc__ or name,'inputSchema':{
        'type':'object','properties':properties,'required':required,'additionalProperties':False}},function


_DEFINITIONS=[
    definition('match',find_matches,['text','left','right'],{
        'text':string(),'left':string(max_length=1000),'right':string(max_length=1000),
        'min_gap':integer(),'max_gap':integer(),'iupac':{'type':'boolean'},
        'strand':string(enum=('+','-','both')),'max_matches':integer(1,1000)}),
    definition('nucleotide_features',nucleotide_features,['text'],{
        'text':string(),'alphabet':string(enum=('dna','rna')),'k':integer(1,6)}),
    definition('rna_fold',fold_rna,['text'],{
        'text':string(max_length=300),'min_loop':integer(0,300),'max_length':integer(1,300),'wobble':{'type':'boolean'}}),
    definition('candidates',generate_candidates,['text'],{
        'text':string(max_length=4096),'mode':string(enum=('preserve','regex','exhaustive','sample','pairwise')),
        'limit':integer(1,1000),'seed':integer(0,2147483647)}),
    definition('compare',compare_aligned,['sequence_a','sequence_b','alphabet'],{
        'sequence_a':string(),'sequence_b':string(),'alphabet':string(enum=('dna','rna','protein')),
        'gap_policy':string(enum=('exclude','mismatch')),'ambiguity_policy':string(enum=('exclude','strict','compatible')),
        'start':integer(1,100000),'end':integer(1,100000)}),
    definition('protein_ionization',ionization,['text'],{
        'text':string(),'ph':{'type':'number','minimum':0,'maximum':14}}),
]
TOOLS={doc['name']:(doc,function) for doc,function in _DEFINITIONS}


def encode(value):
    result=json.dumps(value,ensure_ascii=False,allow_nan=False,separators=(',',':'),sort_keys=True).encode('utf-8')
    if len(result)>MAX_RESPONSE: raise ValueError('response limit exceeded; narrow input/options')
    return result


def call_tool(name,arguments):
    if not isinstance(name,str) or name not in TOOLS: raise ValueError('unknown tool')
    doc,function=TOOLS[name];schema=doc['inputSchema']
    if not isinstance(arguments,dict) or set(arguments)-set(schema['properties']) or not set(schema['required'])<=set(arguments):
        raise ValueError('missing required or unknown argument')
    for key,value in arguments.items():
        spec=schema['properties'][key];kind=spec['type']
        valid=({'string':isinstance(value,str),'integer':type(value) is int,
                'number':type(value) in (int,float),'boolean':type(value) is bool})[kind]
        if not valid: raise ValueError(f'{key}: expected {kind}')
        if isinstance(value,str) and len(value)>spec['maxLength']: raise ValueError(f'{key}: text limit exceeded')
        if 'enum' in spec and value not in spec['enum']: raise ValueError(f'{key}: unsupported value')
        if kind in ('integer','number') and (not spec['minimum']<=value<=spec['maximum'] or not math.isfinite(value)):
            raise ValueError(f'{key}: numeric limit exceeded')
    args=dict(arguments)
    if name=='match':
        args.setdefault('max_matches',1000)
        lo=args.get('min_gap',0);hi=args.get('max_gap',lo)
        if len(args['text'])*max(0,hi-lo+1)*(len(args['left'])+len(args['right']))>10000000:
            raise ValueError('matching work limit exceeded')
    elif name=='candidates':
        args.setdefault('limit',1000)
    result=function(**args)
    encode(result)
    return result


class McpSession:
    def __init__(self): self.state='new'

    def handle(self,request):
        request_id=request.get('id') if isinstance(request,dict) else None
        def error(code,message): return {'jsonrpc':'2.0','id':request_id,'error':{'code':code,'message':message}}
        if (not isinstance(request,dict) or request.get('jsonrpc')!='2.0' or not isinstance(request.get('method'),str)
                or ('id' in request and (type(request_id) not in (int,str)))):
            request_id=None
            return error(-32600,'Invalid Request')
        method=request['method'];params=request.get('params',{})
        if 'id' not in request:
            if method=='notifications/initialized' and self.state=='initializing': self.state='ready'
            return None
        if not isinstance(params,dict): return error(-32602,'params must be an object')
        if method=='initialize':
            if self.state!='new': return error(-32600,'already initialized')
            if not isinstance(params.get('protocolVersion'),str) or not isinstance(params.get('capabilities'),dict) or not isinstance(params.get('clientInfo'),dict):
                return error(-32602,'initialize requires protocolVersion, capabilities and clientInfo')
            self.state='initializing'
            result={'protocolVersion':VERSION,'capabilities':{'tools':{}},'serverInfo':{'name':'sequence-tools','version':'1.0.0'}}
        elif method=='ping': result={}
        elif self.state!='ready': return error(-32600,'initialize and notifications/initialized required')
        elif method=='tools/list':
            if params: return error(-32602,'tools/list accepts no parameters or cursor')
            result={'tools':[doc for doc,_ in TOOLS.values()]}
        elif method=='tools/call':
            if set(params)-{'name','arguments','_meta'} or not isinstance(params.get('name'),str) or params['name'] not in TOOLS:
                return error(-32602,'unknown tool or invalid call parameters')
            try:
                data=call_tool(params['name'],params.get('arguments',{}))
                result={'content':[{'type':'text','text':encode(data).decode('utf-8')}], 'structuredContent':data,'isError':False}
            except (ValueError,TypeError) as failure:
                result={'content':[{'type':'text','text':str(failure)}],'isError':True}
        else: return error(-32601,'Method not found')
        return {'jsonrpc':'2.0','id':request_id,'result':result}


def serve_stdio(input=None,output=None):
    input=sys.stdin.buffer if input is None else input
    output=sys.stdout.buffer if output is None else output
    session=McpSession()
    while True:
        line=input.readline(MAX_REQUEST+1)
        if not line: return
        if len(line)>MAX_REQUEST:
            # Close the session: continuing would interpret a partial line as a new request.
            response={'jsonrpc':'2.0','id':None,'error':{'code':-32700,'message':'request limit exceeded'}}
            output.write(encode(response)+b'\n');output.flush();return
        try:
            request=strict_json(line.decode('utf-8'))
            response=session.handle(request)
        except (ValueError,UnicodeError,RecursionError):
            response={'jsonrpc':'2.0','id':None,'error':{'code':-32700,'message':'Parse error'}}
        if response is not None:
            try: data=encode(response)
            except ValueError: data=encode({'jsonrpc':'2.0','id':response.get('id'),'error':{'code':-32603,'message':'response limit exceeded'}})
            output.write(data+b'\n');output.flush()


class Handler(BaseHTTPRequestHandler):
    def setup(self):
        super().setup()
        self.connection.settimeout(15)

    def do_POST(self):
        if self.path!='/tools/call': self.send_error(404);return
        if self.headers.get('Origin') is not None:
            self.send_error(403,'browser origins are not supported');return
        if self.headers.get_content_type()!='application/json' or self.headers.get('Transfer-Encoding'):
            self.send_error(415);return
        lengths=self.headers.get_all('Content-Length',[])
        try:
            if len(lengths)!=1: raise ValueError('one Content-Length is required')
            length=int(lengths[0])
            if not 0<length<=MAX_REQUEST: raise ValueError('request size limit exceeded')
            raw=self.rfile.read(length)
            if len(raw)!=length: raise ValueError('incomplete request body')
            request=strict_json(raw.decode('utf-8'))
            if not isinstance(request,dict) or set(request)!={'name','arguments'}:
                raise ValueError('body requires name and arguments')
            data=encode(call_tool(request['name'],request['arguments']));status=200
        except (ValueError,UnicodeError,TypeError,RecursionError) as error:
            data=encode({'error':str(error)});status=400
        self.send_response(status);self.send_header('Content-Type','application/json; charset=utf-8')
        self.send_header('Content-Length',str(len(data)));self.end_headers();self.wfile.write(data)

    def log_message(self,*args): pass


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('transport',choices=('stdio','http'))
    parser.add_argument('--host',default='127.0.0.1')
    parser.add_argument('--port',type=int,default=8765)
    args=parser.parse_args(argv)
    if args.transport=='stdio': serve_stdio()
    else:
        with HTTPServer((args.host,args.port),Handler) as server: server.serve_forever()


if __name__=='__main__': main()
