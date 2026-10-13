import gzip
import io
import json
import unittest
from urllib.error import HTTPError
from unittest.mock import patch

from acquisition import Client, FetchError, parse_pubmed, urllib_transport
from sequence_service import McpSession, call_tool, serve_stdio, TOOLS

PUBMED='''<?xml version="1.0"?><!DOCTYPE PubmedArticleSet SYSTEM "https://example.invalid/pubmed.dtd"><PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>123</PMID><Article><ArticleTitle>A <i>nested</i> title</ArticleTitle><Abstract><AbstractText Label="RESULTS">Some <b>evidence</b>.</AbstractText></Abstract></Article></MedlineCitation><PubmedData><ArticleIdList><ArticleId IdType="doi">10.0/example</ArticleId></ArticleIdList></PubmedData></PubmedArticle></PubmedArticleSet>'''


class AcquisitionTests(unittest.TestCase):
    def test_offline_search_pagination_and_cache(self):
        urls=[]
        def transport(url,**kwargs):
            urls.append(url)
            if 'retstart=0' in url: ids=['1','2']
            else: ids=['3']
            return json.dumps({'esearchresult':{'count':'3','idlist':ids}}).encode()
        client=Client(transport=transport,sleep=lambda _:None,cache={})
        pages=list(client.search_pages('a & b',page_size=2,max_records=3))
        self.assertEqual([p['ids'] for p in pages],[['1','2'],['3']])
        self.assertTrue(pages[-1]['complete'])
        self.assertIn('term=a+%26+b',urls[0])
        self.assertEqual(list(client.search_pages('a & b',page_size=2,max_records=3)),pages)
        self.assertEqual(len(urls),2)
        limited=list(client.search_pages('a & b',page_size=2,max_records=2))
        self.assertFalse(limited[-1]['complete'])
        with self.assertRaises(ValueError): client.search('x',start=9999,limit=2)

    def test_retry_limits_and_failure_redaction(self):
        attempts=[];waits=[]
        def transport(url,**kwargs):
            attempts.append(url)
            if len(attempts)<3: raise TimeoutError('secret')
            return b'ok'
        self.assertEqual(Client(transport=transport,sleep=waits.append).fetch('https://example.invalid'),'ok')
        self.assertEqual(waits,[1,2])
        def fail(url,**kwargs): raise HTTPError(url,404,'secret',{},None)
        with self.assertRaisesRegex(FetchError,'HTTPError') as error: Client(transport=fail).fetch('secret')
        self.assertNotIn('secret',str(error.exception))
        with self.assertRaises(FetchError): Client(transport=lambda *a,**k:b'x'*11,max_bytes=10).fetch('url')

    def test_gzip_expansion_cap_and_closed_response(self):
        class Response(io.BytesIO): headers={'Content-Encoding':'gzip'}
        response=Response(gzip.compress(b'A'*5000))
        with patch('acquisition.urlopen',return_value=response):
            with self.assertRaises(FetchError): urllib_transport('https://example.invalid',timeout=1,max_bytes=100)
        self.assertTrue(response.closed)

    def test_pubmed_crossrefs_and_unsupported_xml(self):
        article=parse_pubmed(PUBMED)['articles'][0]
        self.assertEqual(article['title'],'A nested title')
        self.assertEqual(article['abstract'][0]['text'],'Some evidence.')
        self.assertEqual(article['cross_references'][0]['type'],'doi')
        for text in ('<ERROR>bad</ERROR>','<PubmedArticleSet><PubmedBookArticle/></PubmedArticleSet>','<!DOCTYPE x [<!ENTITY a "b">]><x/>'):
            with self.assertRaises(FetchError): parse_pubmed(text)

    def test_fetch_records_are_parsed_and_identifiers_validated(self):
        urls=[]
        def transport(url,**kwargs): urls.append(url);return b'>accession\nACGT\n'
        client=Client(transport=transport,sleep=lambda _:None)
        self.assertEqual(client.sequences(['NM_1.2'])[0].sequence,'ACGT')
        self.assertIn('rettype=fasta',urls[0])
        with self.assertRaises(ValueError): client.sequences(['../secret'])
        with self.assertRaises(ValueError): client.structure('https://example.invalid')


class ServiceTests(unittest.TestCase):
    def request(self,session,method,params=None,id=1):
        return session.handle({'jsonrpc':'2.0','id':id,'method':method,'params':params or {}})

    def initialize(self,session):
        r=self.request(session,'initialize',{'protocolVersion':'2025-11-25','capabilities':{},'clientInfo':{'name':'test','version':'1'}})
        self.assertIn('tools',r['result']['capabilities'])
        self.assertIsNone(session.handle({'jsonrpc':'2.0','method':'notifications/initialized'}))

    def test_handshake_tools_errors_and_notification_silence(self):
        session=McpSession()
        self.assertIn('error',self.request(session,'tools/list'))
        self.initialize(session)
        self.assertEqual(len(self.request(session,'tools/list')['result']['tools']),len(TOOLS))
        result=self.request(session,'tools/call',{'name':'match','arguments':{'text':'ACGT','left':'A','right':'T','min_gap':2}})['result']
        self.assertFalse(result['isError'])
        self.assertEqual(result['structuredContent']['matches'][0]['span'],[1,4])
        bad=self.request(session,'tools/call',{'name':'match','arguments':{'text':'AC!','left':'A','right':'T'}})
        self.assertTrue(bad['result']['isError'])
        self.assertIn('error',self.request(session,'tools/call',{'name':'missing'}))
        self.assertIn('error',self.request(session,'initialize'))
        self.assertIsNone(session.handle({'jsonrpc':'2.0','method':'unknown-notification'}))
        self.assertEqual(session.handle([])['error']['code'],-32600)

    def test_argument_type_work_and_response_limits(self):
        for args in ({'text':'AC','left':'A','right':'C','max_matches':True},
                     {'text':'A'*100000,'left':'A','right':'C','max_gap':10000},
                     {'text':'AC','left':'A','right':'C','filename':'/tmp/secret'}):
            with self.assertRaises(ValueError): call_tool('match',args)
        with self.assertRaises(ValueError): call_tool('protein_ionization',{'text':'AC','ph':float('nan')})
        with self.assertRaises(ValueError): call_tool('rna_fold',{'text':'A'*301})
        with patch('sequence_service.MAX_RESPONSE',20):
            with self.assertRaises(ValueError): call_tool('nucleotide_features',{'text':'ACGT'})

    def test_stdio_bad_json_does_not_stop_following_ping(self):
        source=io.BytesIO(b'not json\n{"jsonrpc":"2.0","id":2,"method":"ping"}\n')
        output=io.BytesIO();serve_stdio(source,output)
        rows=[json.loads(line) for line in output.getvalue().splitlines()]
        self.assertEqual(rows[0]['error']['code'],-32700)
        self.assertEqual(rows[1]['result'],{})
        source=io.BytesIO(b'x'*50+b'\n');output=io.BytesIO()
        with patch('sequence_service.MAX_REQUEST',10): serve_stdio(source,output)
        self.assertEqual(len(output.getvalue().splitlines()),1)

    def test_http_actual_request_boundary(self):
        from http.client import HTTPConnection
        from http.server import HTTPServer
        import threading
        from sequence_service import Handler
        server=HTTPServer(('127.0.0.1',0),Handler)
        thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        try:
            client=HTTPConnection('127.0.0.1',server.server_port,timeout=5)
            body=json.dumps({'name':'match','arguments':{'text':'AT','left':'A','right':'T'}})
            client.request('POST','/tools/call',body,{'Content-Type':'application/json'})
            response=client.getresponse()
            self.assertEqual(response.status,200)
            self.assertEqual(json.loads(response.read())['matches'][0]['span'],[1,2])
            client.close()
            client=HTTPConnection('127.0.0.1',server.server_port,timeout=5)
            client.request('POST','/tools/call',body,{'Content-Type':'application/json','Origin':'https://example.invalid'})
            response=client.getresponse();self.assertEqual(response.status,403);response.read();client.close()
        finally:
            server.shutdown();server.server_close();thread.join(timeout=5)
