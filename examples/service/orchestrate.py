"""Run from repo root: python -S examples/service/orchestrate.py (offline).

A caller may substitute JevClassifier and its own llama.cpp / Pro handlers.
This example does not read keys, call a model, or execute side-effecting tools.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from sequence_router import route_item


def run_pipeline(item, *, classifier=None, threshold=None, local=None, pro=None):
    decision=route_item(item,classifier=classifier,threshold=threshold)
    handler={'local':local,'pro':pro}.get(decision['route'])
    if handler is None:
        return {'decision':decision,'result':None}
    try:
        return {'decision':decision,'result':handler(item)}
    except Exception as error:
        return {'decision':decision,'result':None,'fallback':'human','error_type':type(error).__name__}


if __name__=='__main__':
    import json
    item={'kind':'literature','objective':'inspect evidence','title':'Example title'}
    print(json.dumps(run_pipeline(item,local=lambda item:{'queued_for_local_analysis':True}),sort_keys=True))
