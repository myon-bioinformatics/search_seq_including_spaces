"""Run from repo root: python -S -m scripts.benchmark_iupac_masks.

Load a fixed pre-change reference using git; compare complete outputs before timing.
"""
import json, platform, random, statistics, subprocess, sys, time, types

import anchor_matcher
import streaming_matcher

BASELINE = '5ab1885919c52c46a255c7944ef475ae8061c46d'


def reference(name, anchor_reference=None):
    module = types.ModuleType('_baseline_' + name)
    module.__file__ = name + '.py'
    sys.modules[module.__name__] = module
    code = subprocess.check_output(['git', 'show', BASELINE + ':' + name + '.py'])
    # Streaming imports must resolve to the reference anchor module.
    original = sys.modules['anchor_matcher']
    if name == 'streaming_matcher':
        sys.modules['anchor_matcher'] = anchor_reference
    try:
        exec(compile(code, module.__file__, 'exec'), module.__dict__)
    finally:
        sys.modules['anchor_matcher'] = original
    return module



def main():
    old_anchor = reference("anchor_matcher")
    old_stream = reference("streaming_matcher", old_anchor)
    rng = random.Random(20261104)
    text = ''.join(rng.choice('ACGT') for _ in range(20000))
    options = dict(iupac=True, strand='both', max_gap=2, max_matches=10000)
    rows = []
    for name, before, after in (
        ('whole', lambda: old_anchor.find_matches(text, 'ACR', 'YGT', **options),
         lambda: anchor_matcher.find_matches(text, 'ACR', 'YGT', **options)),
        ('chunks_1024', lambda: list(old_stream.match_chunks(
            (text[i:i+1024] for i in range(0, len(text), 1024)), 'ACR', 'YGT', **options)),
         lambda: list(streaming_matcher.match_chunks(
            (text[i:i+1024] for i in range(0, len(text), 1024)), 'ACR', 'YGT', **options))),
    ):
        if before() != after():
            raise SystemExit(f"output mismatch: {name}")
        times = [[], []]
        for turn in range(9):
            for index in ((0, 1) if turn % 2 == 0 else (1, 0)):
                start = time.perf_counter_ns()
                (before, after)[index]()
                times[index].append(time.perf_counter_ns() - start)
        rows.append(dict(case=name, before_ns=times[0], after_ns=times[1],
                         before_median_ns=statistics.median(times[0]),
                         after_median_ns=statistics.median(times[1]), output_equal=True))
    print(json.dumps(dict(baseline=BASELINE, python=platform.python_version(),
                         platform=platform.platform(), seed=20261104, sequence_length=len(text),
                         repeats=9, measurements=rows), indent=2))


if __name__ == '__main__':
    main()
