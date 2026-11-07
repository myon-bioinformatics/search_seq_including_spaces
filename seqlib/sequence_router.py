"""Optional injected decision boundary for Jev/local/pro routing, without API I/O."""
import math
from collections.abc import Mapping

QUESTIONS = {'sequence_type': ('dna','rna','protein','unknown'),
             'is_relevant_species': ('yes','no','unknown'),
             'is_pseudogene_candidate': ('yes','no','unknown'),
             'needs_deeper_analysis': ('yes','no','unknown')}


def route_item(item, *, question='needs_deeper_analysis', classifier=None,
               threshold=None, fallback='local', important=False):
    """classifier(request) returns {decision,confidence}; any failure falls back.

    Confidence is a routing score, never a calibrated probability. No action,
    network call or generation is executed here. Callers own transport timeout.
    """
    if not isinstance(item,Mapping) or item.get('kind') not in ('sequence','literature','analysis'):
        raise ValueError('item kind must be sequence, literature or analysis')
    if question not in QUESTIONS or fallback not in ('local','pro','human') or not isinstance(important,bool):
        raise ValueError('invalid question/fallback/important option')
    if classifier is not None and (isinstance(threshold,bool) or not isinstance(threshold,(int,float)) or not math.isfinite(threshold) or not 0<=threshold<=1):
        raise ValueError('injected classifier requires an explicit evaluated threshold in [0,1]')
    result={'schema':'sequence-route/1','kind':item['kind'],'question':question,'decision':'unknown',
            'confidence':None,'threshold':threshold,'route':'pro' if important else fallback,
            'reason':'important' if important else 'no-classifier','classifier_used':False}
    if important: return result
    if question=='sequence_type' and item['kind']=='sequence':
        raw=item.get('sequence')
        if not isinstance(raw,str): raise ValueError('sequence item requires text')
        sequence=''.join(raw.split()).upper()
        # A/C/G-only strings fit both nucleotide alphabets: leave uncertainty explicit.
        if item.get('alphabet')=='dna' and sequence and set(sequence)<=set('ACGTRYSWKMBDHVN'):
            result.update(decision='dna',route='deterministic',reason='alphabet-rule');return result
        if item.get('alphabet')=='rna' and sequence and set(sequence)<=set('ACGURYSWKMBDHVN'):
            result.update(decision='rna',route='deterministic',reason='alphabet-rule');return result
    if classifier is None: return result
    try:
        response=classifier({'question':question,'choices':list(QUESTIONS[question]),'item':dict(item)})
        if not isinstance(response,dict) or set(response)!={'decision','confidence'} or response['decision'] not in QUESTIONS[question]:
            raise ValueError('invalid decision schema')
        confidence=response['confidence']
        if isinstance(confidence,bool) or not isinstance(confidence,(int,float)) or not math.isfinite(confidence) or not 0<=confidence<=1:
            raise ValueError('invalid confidence')
        result.update(decision=response['decision'],confidence=confidence,classifier_used=True,reason='uncertain')
        if confidence>=threshold and response['decision']!='unknown':
            # Exclusion is only a returned routing recommendation, never an action.
            result.update(route='exclude' if response['decision']=='no' else 'local',reason='threshold-met')
    except Exception as error:
        result.update(classifier_used=True,reason='classifier-failed',error_type=type(error).__name__)
    return result
