"""Optional stdlib acquisition. No network access during import or core analysis."""
import gzip
import hashlib
import io
import math
import re
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

from result_summary import integer
from sequence_readers import read_sequences
from structured_readers import strict_json


class FetchError(ValueError):
    pass


def urllib_transport(url, *, timeout, max_bytes):
    request = Request(url, headers={'User-Agent': 'search-seq-including-spaces/1', 'Accept-Encoding': 'gzip'})
    with urlopen(request, timeout=timeout) as response:
        body = response.read(max_bytes + 1)
        if len(body) > max_bytes:
            raise FetchError('compressed/received response exceeds byte limit')
        if response.headers.get('Content-Encoding') == 'gzip' or body.startswith(b'\x1f\x8b'):
            with gzip.GzipFile(fileobj=io.BytesIO(body)) as compressed:
                body = compressed.read(max_bytes + 1)
        if len(body) > max_bytes:
            raise FetchError('expanded response exceeds byte limit')
        return body


class Client:
    """Injected transport(url, timeout=..., max_bytes=...) -> decoded bytes.

    Optional dict-like cache stores exact response bytes, keyed by URL hash.
    NCBI requests are paced at <=3/s even with an API key. Retries only cover
    transient transport failures and 429/5xx; sizes and parsing fail immediately.
    """
    def __init__(self, *, transport=urllib_transport, timeout=15, max_bytes=8*1024*1024,
                 retries=2, sleep=time.sleep, clock=time.monotonic, cache=None, tool='search_seq_including_spaces', email=None):
        if isinstance(timeout, bool) or not isinstance(timeout, (int,float)) or not math.isfinite(timeout) or timeout <= 0:
            raise ValueError('timeout must be finite and positive')
        integer(max_bytes,'max_bytes',1)
        integer(retries,'retries')
        if retries > 5 or not re.fullmatch(r'[A-Za-z0-9_.-]+',tool):
            raise ValueError('retries must be <= 5 and tool must be a simple identifier')
        self.transport,self.timeout,self.max_bytes,self.retries=transport,timeout,max_bytes,retries
        self.sleep,self.clock,self.cache=sleep,clock,cache
        self.tool,self.email,self.last_request=tool,email,None

    def fetch(self,url,*,pace=False):
        key=hashlib.sha256(url.encode()).hexdigest()
        if self.cache is not None and key in self.cache:
            body=self.cache[key]
        else:
            for attempt in range(self.retries+1):
                if pace and self.last_request is not None:
                    self.sleep(max(0, 0.34-(self.clock()-self.last_request)))
                self.last_request=self.clock()
                try:
                    body=self.transport(url,timeout=self.timeout,max_bytes=self.max_bytes)
                    break
                except (HTTPError,URLError,TimeoutError,OSError) as error:
                    retry = not isinstance(error,HTTPError) or error.code==429 or 500<=error.code<600
                    if not retry or attempt==self.retries:
                        if isinstance(error,HTTPError): error.close()
                        # Do not echo query URLs, email addresses or transport-specific secrets.
                        raise FetchError(f'fetch failed ({type(error).__name__})') from None
                    delay=2**attempt
                    if isinstance(error,HTTPError):
                        retry_after=error.headers.get('Retry-After','') if error.headers else ''
                        if retry_after.isdigit(): delay=min(60,max(delay,int(retry_after)))
                        error.close()
                    self.sleep(delay)
            if not isinstance(body,bytes) or len(body)>self.max_bytes:
                raise FetchError('transport returned invalid or oversized bytes')
            if self.cache is not None:
                self.cache[key]=body
        if not isinstance(body,bytes) or len(body)>self.max_bytes:
            raise FetchError('cache returned invalid or oversized bytes')
        return body.decode('utf-8')

    def _ncbi(self,utility,**params):
        params['tool']=self.tool
        if self.email: params['email']=self.email
        return self.fetch('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'+utility+'.fcgi?'+urlencode(params),pace=True)

    @staticmethod
    def _database(db):
        if db not in ('pubmed','nuccore','protein','structure'):
            raise ValueError('unsupported NCBI database')

    def search(self,term,*,db='pubmed',start=0,limit=100):
        self._database(db)
        integer(start,'start');integer(limit,'limit',1)
        if limit>1000 or not isinstance(term,str) or not term.strip():
            raise ValueError('search needs nonempty term and limit <= 1000')
        if db=='pubmed' and start+limit>10000:
            raise ValueError('PubMed ESearch window exceeds 10000; partition the query')
        data=strict_json(self._ncbi('esearch',db=db,term=term,retstart=start,retmax=limit,retmode='json'))
        result=data.get('esearchresult',{})
        if 'error' in data or result.get('errorlist') or 'count' not in result or not isinstance(result.get('idlist'),list):
            raise FetchError('invalid/failed ESearch response')
        count=int(result['count']);ids=result['idlist']
        if count<0 or len(ids)>limit or any(not isinstance(i,str) or not i.isdigit() for i in ids):
            raise FetchError('invalid ESearch count/ids')
        return {'schema':'ncbi-search/1','database':db,'ids':ids,'count':count,'start':start,
                'next_start':start+len(ids) if start+len(ids)<count else None,
                'complete':start+len(ids)>=count}

    def search_pages(self,term,*,db='pubmed',page_size=100,max_records=1000):
        integer(max_records,'max_records',1)
        integer(page_size,'page_size',1)
        start=0
        seen=set()
        while start<max_records:
            page=self.search(term,db=db,start=start,limit=min(page_size,max_records-start))
            if not page['ids'] and not page['complete']:
                raise FetchError('ESearch made no pagination progress')
            if seen.intersection(page['ids']) or len(set(page['ids']))!=len(page['ids']):
                raise FetchError('ESearch repeated ids across pages; retry with a stable query')
            seen.update(page['ids'])
            yield page
            if page['complete']: return
            start=page['next_start']
        # The last page keeps complete=false / next_start for an explicit resumable bound.

    @staticmethod
    def _ids(ids):
        ids=list(ids)
        if not 1<=len(ids)<=200 or any(not isinstance(i,str) or not re.fullmatch(r'[A-Za-z0-9_.]+',i) for i in ids):
            raise ValueError('provide 1..200 simple accession/UID strings')
        return ','.join(ids)

    def summary(self,ids,*,db='pubmed'):
        self._database(db)
        data=strict_json(self._ncbi('esummary',db=db,id=self._ids(ids),retmode='json'))
        if 'error' in data or not isinstance(data.get('result'),dict):
            raise FetchError('invalid/failed ESummary response')
        return {'schema':'ncbi-summary/1','database':db,'result':data['result']}

    def sequences(self,ids,*,db='nuccore',format='fasta'):
        if db not in ('nuccore','protein') or format not in ('fasta','genbank'):
            raise ValueError('sequences supports nuccore/protein and fasta/genbank')
        rettype='fasta' if format=='fasta' else 'gb' if db=='nuccore' else 'gp'
        text=self._ncbi('efetch',db=db,id=self._ids(ids),rettype=rettype,retmode='text')
        return list(read_sequences(io.StringIO(text),format=format,source=f'NCBI:{db}'))

    def pubmed(self,ids):
        text=self._ncbi('efetch',db='pubmed',id=self._ids(ids),retmode='xml')
        return parse_pubmed(text)

    def structure(self,identifier,*,format='cif'):
        if not isinstance(identifier,str) or not re.fullmatch(r'[0-9][A-Za-z0-9]{3}',identifier) or format not in ('cif','pdb'):
            raise ValueError('RCSB subset requires a legacy four-character PDB id and cif/pdb format')
        url=f'https://files.rcsb.org/download/{identifier.upper()}.{format}'
        text=self.fetch(url)
        return {'schema':'structure-download/1','id':identifier.upper(),'format':format,'text':text,
                'source':url,'sha256':hashlib.sha256(text.encode()).hexdigest(),
                'mmdb_lookup':'https://www.ncbi.nlm.nih.gov/structure/?term='+identifier.upper()}


def parse_pubmed(text):
    # NCBI's external DTD declaration is allowed; internal entity declarations are not.
    if '<!ENTITY' in text.upper() or re.search(r'<!DOCTYPE[^>]*\[',text,re.I):
        raise FetchError('internal XML declarations are unsupported')
    try:
        root=ET.fromstring(text)
    except ET.ParseError as error:
        raise FetchError(f'invalid PubMed XML: {error}') from None
    if root.tag!='PubmedArticleSet' or root.find('.//ERROR') is not None:
        raise FetchError('expected PubmedArticleSet')
    result=[]
    def content(node): return ''.join(node.itertext()) if node is not None else ''
    for entry in root:
        if entry.tag!='PubmedArticle':
            raise FetchError('PubMed parser subset supports PubmedArticle, not book records')
        citation=entry.find('MedlineCitation')
        if citation is None or not citation.findtext('PMID'):
            raise FetchError('article missing PMID')
        result.append({'pmid':citation.findtext('PMID'),
                       'title':content(citation.find('Article/ArticleTitle')),
                       'abstract':[{'label':node.get('Label'),'text':content(node)} for node in citation.findall('Article/Abstract/AbstractText')],
                       'cross_references':[{'type':node.get('IdType'),'id':content(node)} for node in entry.findall('PubmedData/ArticleIdList/ArticleId')]})
    return {'schema':'pubmed-records/1','articles':result}
