#!/usr/bin/env python3
"""Figure 5: the Aldous-Shepp sandwich 1/n <= CV^2_min(n,R) <= 1/(n-1).
Lower bound = Erlang (Aldous-Shepp); upper bound = one-trap construction; the
star is loaded from and audited against data/optimum_n3_R2.json. Reproduces
fig5_sandwich.png and fig5_sandwich.pdf.
"""
import json
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from tempcomp import fpt
from tempcomp.oscillators import one_trap_construction_subgen
from tempcomp.optimum import audit_candidate, load_candidate

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT_PNG = os.path.join(ROOT, "figures", "fig5_sandwich.png")
OUT_PDF = os.path.join(ROOT, "figures", "fig5_sandwich.pdf")
OPT = os.path.join(ROOT, "data", "optimum_n3_R2.json")


def main():
    # exact construction points
    cn, cv = [], []
    for n in [3, 4, 5, 6]:
        cn.append(n)
        cv.append(fpt.cv2(one_trap_construction_subgen(1.0, n)))
    print("construction CV^2:", {n: round(c, 4) for n, c in zip(cn, cv)})

    opt = load_candidate(OPT)
    audit = audit_candidate(opt)
    print("audited cached candidate:", json.dumps(audit, sort_keys=True))

    ns = np.arange(2, 11)
    plt.figure(figsize=(7.6, 5))
    plt.plot(ns, 1 / ns, "o-", color="#16a085", lw=2,
             label=r"$1/n$  (Aldous--Shepp floor)")
    plt.plot(ns, 1 / (ns - 1), "s-", color="#c0392b", lw=2,
             label=r"$1/(n-1)$  (one-trap construction)")
    plt.fill_between(ns, 1 / ns, 1 / (ns - 1), color="orange", alpha=0.15)
    plt.plot(cn, cv, "D", color="#c0392b", ms=7, mfc="white",
             label="construction (exact)")
    plt.plot([opt["n"]], [audit["CV2"]], "*", color="k", ms=16,
             label="audited numerical candidate")
    plt.xlabel("number of internal states  n")
    plt.ylabel(r"min CV$^2$ of a compensating dwell")
    plt.legend(fontsize=9, loc="upper right")
    plt.ylim(0, 0.6)
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=300)
    plt.savefig(OUT_PDF)
    print("wrote", OUT_PNG, "and", OUT_PDF)


if __name__ == "__main__":
    main()
