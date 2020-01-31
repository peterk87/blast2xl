from typing import Mapping, Dict, Any, Set


def invert_dict(d: Mapping) -> Dict[Any, Set]:
    out = {}
    for k, v in d.items():
        if v in out:
            out[v].add(k)
        else:
            out[v] = {k}
    return out
