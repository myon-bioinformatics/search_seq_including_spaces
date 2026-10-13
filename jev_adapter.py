"""Opt-in Jev Choice adapter; official API shape, no SDK dependency.

No key lookup or I/O at import. Pass an injected transport for offline tests or
an explicit key for live requests. The router retains final routing authority.
"""
import json
import math
from urllib.error import URLError, HTTPError
from urllib.request import Request, build_opener, HTTPRedirectHandler

from result_summary import integer
from sequence_router import QUESTIONS
from structured_readers import strict_json

ENDPOINT = 'https://api.typesafe.ai/v1/systemone'
INSTRUCTIONS = {
    'sequence_type':'Select the sequence alphabet supported by the provided evidence. Return unknown when multiple alphabets fit.',
    'is_relevant_species':'Does the supplied evidence match the explicitly requested species? Return unknown without sufficient metadata.',
    'is_pseudogene_candidate':'Do the supplied deterministic features justify follow-up as a pseudogene candidate? This is screening, not a scientific conclusion.',
    'needs_deeper_analysis':'Does the supplied sequence, literature or analysis need further analysis for the stated objective? Return unknown if the objective or evidence is insufficient.'}


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self,req,fp,code,msg,headers,newurl):
        raise URLError('redirect refused')


class JevClassifier:
    def __init__(self, *, model, transport=None, api_key=None, timeout=15, max_bytes=1024*1024):
        if not isinstance(model,str) or not model.strip(): raise ValueError('explicit model is required')
        if isinstance(timeout,bool) or not isinstance(timeout,(int,float)) or not math.isfinite(timeout) or timeout<=0:
            raise ValueError('timeout must be finite and positive')
        integer(max_bytes,'max_bytes',1)
        if transport is None and (not isinstance(api_key,str) or not api_key or '\n' in api_key or '\r' in api_key):
            raise ValueError('live Jev requires an explicitly supplied API key')
        self.model,self.timeout,self.max_bytes=model,timeout,max_bytes
        self.transport=transport
        self._api_key=api_key

    def _request(self,payload):
        data=json.dumps(payload,ensure_ascii=False,allow_nan=False).encode('utf-8')
        if len(data)>self.max_bytes: raise ValueError('Jev request size limit exceeded')
        if self.transport is not None:
            body=self.transport(payload,timeout=self.timeout,max_bytes=self.max_bytes)
        else:
            request=Request(ENDPOINT,data=data,headers={'Authorization':'Bearer '+self._api_key,'Content-Type':'application/json'})
            try:
                with build_opener(_NoRedirect()).open(request,timeout=self.timeout) as response:
                    body=response.read(self.max_bytes+1)
            except (OSError,HTTPError) as error:
                if isinstance(error,HTTPError): error.close()
                raise ValueError(f'Jev request failed ({type(error).__name__})') from None
        if not isinstance(body,bytes) or len(body)>self.max_bytes: raise ValueError('Jev response size/type limit exceeded')
        return strict_json(body.decode('utf-8'))

    def __call__(self,request):
        question=request.get('question')
        if question not in QUESTIONS or request.get('choices')!=list(QUESTIONS[question]):
            raise ValueError('unsupported atomic routing question or choices')
        payload={'model':self.model,'state':request['item'],'questions':{
            question:{'type':'choice','instructions':INSTRUCTIONS[question],
                      'criteria':{choice:None for choice in QUESTIONS[question]}}}}
        response=self._request(payload)
        try:
            answer=response['answers'][question]
            probabilities=answer['probabilities']
            if answer['type']!='choice' or set(probabilities)!=set(QUESTIONS[question]): raise ValueError('unexpected Choice fields')
            if any(type(p) not in (int,float) or not math.isfinite(p) or not 0<=p<=1 for p in probabilities.values()):
                raise ValueError('invalid Choice probabilities')
            if not math.isclose(sum(probabilities.values()),1.0,abs_tol=1e-6): raise ValueError('Choice probabilities do not sum to one')
            choice,confidence=answer['choice'],answer['confidence']
            if (choice not in probabilities or probabilities[choice]!=max(probabilities.values()) or
                    type(confidence) not in (int,float) or not math.isfinite(confidence) or not 0<=confidence<=1):
                raise ValueError('invalid Choice decision/confidence')
        except (KeyError,TypeError,AttributeError) as error:
            raise ValueError('invalid Jev response schema') from None
        return {'decision':choice,'confidence':confidence}
