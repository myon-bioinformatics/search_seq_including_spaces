"""Routing boundary (seqlib/sequence_router.py) and its optional Jev classifier (seqlib/jev_adapter.py)."""
import json, unittest
from seqlib.jev_adapter import JevClassifier
from seqlib.sequence_router import route_item


class RouterTests(unittest.TestCase):
    def test_routing_failure_uncertainty_no_side_effect(self):
        item={'kind':'analysis','features':{'score':4}}
        self.assertEqual(route_item(item)['route'],'local')
        def fail(request): raise TimeoutError('private transport details')
        r=route_item(item,classifier=fail,threshold=0.8)
        self.assertEqual(r['reason'],'classifier-failed')
        self.assertNotIn('private',str(r))
        good=lambda request:{'decision':'no','confidence':0.9}
        self.assertEqual(route_item(item,classifier=good,threshold=0.8)['route'],'exclude')
        self.assertEqual(route_item(item,classifier=good,threshold=0.95)['route'],'local')
        self.assertEqual(route_item(item,classifier=fail,threshold=0.8,important=True)['route'],'pro')
        with self.assertRaises(ValueError): route_item(item,classifier=good)
        r=route_item({'kind':'sequence','sequence':'ACGT','alphabet':'dna'},question='sequence_type',classifier=fail,threshold=0.8)
        self.assertEqual(r['route'],'deterministic')


class JevTests(unittest.TestCase):
    def test_official_choice_wire_contract_and_routing(self):
        requests=[]
        def transport(payload,**options):
            requests.append((payload,options))
            return json.dumps({'answers':{'needs_deeper_analysis':{'type':'choice','choice':'yes',
                'confidence':0.8,'probabilities':{'yes':0.9,'no':0.05,'unknown':0.05}}}}).encode()
        classifier=JevClassifier(model='test-model',transport=transport,timeout=3)
        result=route_item({'kind':'analysis','data':'fixture'},classifier=classifier,threshold=0.7)
        self.assertEqual(result['route'],'local')
        self.assertEqual(result['reason'],'threshold-met')
        self.assertEqual(requests[0][0]['questions']['needs_deeper_analysis']['type'],'choice')
        self.assertEqual(requests[0][1]['timeout'],3)

    def test_invalid_response_timeout_and_size_fall_back(self):
        for body in (b'{}',b'not json',b'x'*101):
            classifier=JevClassifier(model='test',transport=lambda *a,**k:body,max_bytes=1000)
            self.assertEqual(route_item({'kind':'analysis'},classifier=classifier,threshold=0.8)['reason'],'classifier-failed')
        def timeout(*args,**kwargs): raise TimeoutError('private details')
        classifier=JevClassifier(model='test',transport=timeout)
        result=route_item({'kind':'literature'},classifier=classifier,threshold=0.8)
        self.assertEqual(result['route'],'local')
        self.assertNotIn('private',str(result))
        with self.assertRaises(ValueError): JevClassifier(model='test')
