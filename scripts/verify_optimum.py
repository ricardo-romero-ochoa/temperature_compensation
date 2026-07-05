#!/usr/bin/env python3
"""Fast audit of the cached Figure 5 constrained-optimization candidate."""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from tempcomp.optimum import audit_candidate, load_candidate

DEFAULT = os.path.join(os.path.dirname(__file__), "..", "data", "optimum_n3_R2.json")


def main() -> None:
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    candidate = load_candidate(path)
    audit = audit_candidate(candidate)
    print(json.dumps(audit, indent=2, sort_keys=True))
    print("cached optimum candidate verified")


if __name__ == "__main__":
    main()
