#!/usr/bin/env python3
"""Compare two SLUCM single-column trajectories.

Item 6 of dev/plan-slucm-coupling.md section 8 sets the gate at 1e-10 relative
on SH/LH/G/TS/TR/TB/TG. That is the *allowance*, not the target: the two builds
compile the same Fortran with the same promotion flags, so the honest
expectation is bit-identical output and any non-zero difference is a finding.
Both numbers are therefore reported -- whether the trajectories are bitwise
equal, and, if not, the worst relative difference against the gate.

    ./compare_trajectories.py trajectory_erf.csv trajectory_wrf.csv

Exit status 0 only when every gated column passes.
"""

import csv
import math
import sys

GATE = 1.0e-10
GATED = ["SH", "LH", "G", "TS", "TR", "TB", "TG"]


def load(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    if not rows:
        raise SystemExit(f"{path}: empty")
    head = [h.strip() for h in rows[0]]
    data = [[float(x) for x in r] for r in rows[1:] if r and r[0].strip()]
    return head, data


def main(argv):
    if len(argv) != 3:
        raise SystemExit(__doc__)
    ha, da = load(argv[1])
    hb, db = load(argv[2])

    if ha != hb:
        print("FAIL: column headers differ")
        return 1
    if len(da) != len(db):
        print(f"FAIL: step count differs -- {len(da)} vs {len(db)}")
        return 1
    if not da:
        print("FAIL: no rows to compare")
        return 1

    print(f"steps compared : {len(da)}")
    print(f"columns        : {len(ha)}")
    print()

    bitwise = True
    failed = []
    widest = 0.0

    for j, name in enumerate(ha):
        if name in ("step", "time_h"):
            continue
        worst = 0.0
        worst_step = -1
        exact = True
        for i, (ra, rb) in enumerate(zip(da, db)):
            x, y = ra[j], rb[j]
            if x != y:
                exact = False
                denom = max(abs(x), abs(y))
                rel = 0.0 if denom == 0.0 else abs(x - y) / denom
                if rel > worst:
                    worst, worst_step = rel, i + 1
        if not exact:
            bitwise = False
        gated = name in GATED
        if gated:
            widest = max(widest, worst)
            if worst > GATE:
                failed.append((name, worst, worst_step))
        mark = "identical" if exact else f"max rel {worst:.3e} @ step {worst_step}"
        flag = " [gated]" if gated else ""
        print(f"  {name:<6} {mark}{flag}")

    print()
    if bitwise:
        print("RESULT: BITWISE_IDENTICAL")
        print("The vendored copy and upstream produce identical bits in this")
        print("configuration -- the ERF-local edits are numerically inert here.")
        return 0

    if failed:
        print(f"RESULT: FAIL -- {len(failed)} gated column(s) exceed {GATE:.0e}")
        for n, w, s in failed:
            print(f"  {n}: {w:.3e} at step {s}")
        return 1

    print(f"RESULT: WITHIN_GATE (worst gated {widest:.3e} < {GATE:.0e})")
    print("Not bitwise. Under identical sources and flags that is itself worth")
    print("explaining before this is accepted as a pass.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
