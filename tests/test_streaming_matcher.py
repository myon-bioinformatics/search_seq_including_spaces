import io
import itertools
import random
import unittest
from anchor_matcher import find_matches
from streaming_matcher import match_chunks, match_large_file


def key(match):
    return (match['strand'], match['span'], match['gap'])


class ChunkTests(unittest.TestCase):
    def test_randomized_oracle_across_all_chunk_boundaries(self):
        rng=random.Random(13)
        for _ in range(35):
            text=''.join(rng.choice('ACGTRN \n') for _ in range(40))
            for strand,iupac in itertools.product(('+','-','both'),(True,False)):
                options=dict(min_gap=0,max_gap=3,strand=strand,iupac=iupac,max_matches=10000)
                expected=sorted(find_matches(text,'A','T',**options)['matches'],key=key)
                for size in (1,3,9,64):
                    events=list(match_chunks((text[i:i+size] for i in range(0,len(text),size)),'A','T',**options))
                    self.assertEqual(sorted([e['match'] for e in events[:-1]],key=key),expected)
                    self.assertEqual(events[-1]['match_count'],len(expected))
                    self.assertEqual(events[-1]['source_length'],len(text))

    def test_multibase_anchors_spanning_chunks(self):
        text=' AA C\nNN GTT AACNNGTT'
        for size in range(1,12):
            chunks=(text[i:i+size] for i in range(0,len(text),size))
            report=list(match_chunks(chunks,'AAC','GTT',min_gap=2,strand='both'))
            expected=find_matches(text,'AAC','GTT',min_gap=2,strand='both')
            self.assertEqual(sorted([r['match'] for r in report[:-1]],key=key),sorted(expected['matches'],key=key))

    def test_limits_empty_invalid_input_never_emit_completion(self):
        events=match_chunks(['AAAA'],'A','A',max_matches=1)
        self.assertEqual(next(events)['event'],'match')
        with self.assertRaises(ValueError): list(events)
        with self.assertRaises(ValueError): list(match_chunks(['AC!'],'A','C'))
        with self.assertRaises(ValueError): list(match_chunks([],'A','C',max_gap=100,max_span=10))
        self.assertEqual(list(match_chunks([],'A','C'))[-1]['match_count'],0)

    def test_fasta_record_boundaries_comments_and_handle_ownership(self):
        handle=io.StringIO('>x\nA\n;ignored\nC\n>y\nAC\n')
        events=list(match_large_file(handle,'A','C',chunk_size=8))
        self.assertEqual([(e['record']['id'],e['event']) for e in events], [('x','match'),('x','complete'),('y','match'),('y','complete')])
        self.assertFalse(handle.closed)
        self.assertEqual(events[0]['match']['source_span'],[1,3])
        handle=io.StringIO('>x\nA\n>y\nC\n')
        self.assertEqual([e['match_count'] for e in match_large_file(handle,'A','C',chunk_size=8)],[0,0])

    def test_long_unwrapped_lines_use_bounded_reads(self):
        class Guard(io.StringIO):
            def readline(self,size=-1):
                if size<0 or size>17: raise AssertionError('unbounded read')
                return super().readline(size)
        handle=Guard('>x\n'+'C'*10000+'AT\n')
        events=list(match_large_file(handle,'A','T',chunk_size=16))
        self.assertEqual(events[0]['match']['span'],[10001,10002])
        self.assertEqual(events[-1]['sequence_length'],10002)


    def test_fastq_chunks_quality_and_completion(self):
        for size in (4,8,17):
            events=list(match_large_file(io.StringIO('@x\nACGTACGT\n+\n!!!!!!!!\n@y\nACGT\n+\n@@@@\n'),'A','T',min_gap=2,format='fastq',chunk_size=size))
            self.assertEqual([e['match_count'] for e in events if e['event']=='complete'],[2,1])
        for text in ('@x\nACGT\n+\n!!!\n','@x\nACGT\n+\n!!!!!\n','@x\nAC GT\n+\n!!!!!\n'):
            events=match_large_file(io.StringIO(text),'A','T',min_gap=2,format='fastq',chunk_size=8)
            received=[]
            with self.assertRaises(ValueError):
                for event in events: received.append(event)
            self.assertFalse(any(e['event']=='complete' for e in received))

    def test_fastq_crlf_split_chunks(self):
        for size in (5,6,7,8):
            text='@x\r\nACGTACGT\r\n+\r\n!!!!!!!!\r\n'
            events=list(match_large_file(io.StringIO(text),'A','T',min_gap=2,format='fastq',chunk_size=size))
            self.assertEqual(events[-1]['match_count'],2)
