#!/usr/bin/env python3
"""Compare two SLUCM smoke work directories field by field.

The README asserts that one rank and two ranks agree *bit for bit* on every
Noah-MP land field, and calls that "the decomposition-independence check worth
having". But ``run_smoke.sh`` runs a single decomposition per invocation and
``check_smoke.py`` only ever looks inside one work directory, so nothing in the
tree actually checked it -- the claim rested on a comparison someone did by
hand once. This closes that gap:

    ./run_smoke.sh --exe ... --work /tmp/np1 --np 1
    ./run_smoke.sh --exe ... --work /tmp/np2 --np 2
    ./compare_ranks.py /tmp/np1 /tmp/np2

Exit status is 0 only when every field in every case is identical.

Bit-for-bit is the right bar here and a deliberately harsh one: a domain
decomposition must not change a per-column land surface calculation at all. A
difference of any size is a bug -- in halo exchange, in the Noah-MP box
assumptions, or in the urban tile's indexing -- and not something to be waved
through with a tolerance.
"""

import argparse
import glob
import os
import sys

import numpy as np
from netCDF4 import Dataset

# run_smoke.sh names the two configurations this way: `rural` is the
# bulk-urban control, not an absence of urban.
CASES = ("urban", "rural")


def land_file(work, case):
    """Latest lnd<step>/Level_0.nc for one case, or None if the case is absent."""
    hits = sorted(glob.glob(os.path.join(work, case, "lnd*", "Level_0.nc")))
    return hits[-1] if hits else None


def identical(a, b):
    """Exact equality, treating NaN as equal to NaN.

    Plain ``==`` reports NaN != NaN, which would flag every masked land point
    as a difference. Land output is full of them, so that distinction matters.
    """
    if a.shape != b.shape:
        return False, f"shape {a.shape} vs {b.shape}"
    if a.dtype.kind == "f" or b.dtype.kind == "f":
        af = np.asarray(a, dtype="float64")
        bf = np.asarray(b, dtype="float64")
        same = (af == bf) | (np.isnan(af) & np.isnan(bf))
        if same.all():
            return True, ""
        d = np.abs(af - bf)
        d = d[np.isfinite(d)]
        worst = d.max() if d.size else float("nan")
        return False, f"{(~same).sum()} of {same.size} points differ, max|diff|={worst:.6e}"
    if np.array_equal(a, b):
        return True, ""
    return False, f"{(a != b).sum()} of {a.size} points differ"


def compare_case(fa, fb):
    """Return (n_compared, [(field, reason), ...])."""
    diffs = []
    with Dataset(fa) as da, Dataset(fb) as db:
        names = sorted(set(da.variables) | set(db.variables))
        for n in names:
            if n not in da.variables:
                diffs.append((n, "present only in the second run"))
                continue
            if n not in db.variables:
                diffs.append((n, "present only in the first run"))
                continue
            ok, why = identical(np.asarray(da[n][...]), np.asarray(db[n][...]))
            if not ok:
                diffs.append((n, why))
    return len(names), diffs


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("work_a", help="work directory of the first run (e.g. --np 1)")
    p.add_argument("work_b", help="work directory of the second run (e.g. --np 2)")
    args = p.parse_args(argv)

    failed = False
    checked_any = False

    for case in CASES:
        fa = land_file(args.work_a, case)
        fb = land_file(args.work_b, case)
        if fa is None or fb is None:
            missing = args.work_a if fa is None else args.work_b
            print(f"FAIL {case}: no lnd*/Level_0.nc under {missing}")
            failed = True
            continue

        checked_any = True
        n, diffs = compare_case(fa, fb)
        if diffs:
            failed = True
            print(f"FAIL {case}: {len(diffs)} of {n} fields differ")
            for name, why in diffs:
                print(f"       {name}: {why}")
        else:
            print(f"ok   {case}: all {n} fields bit-identical")

    if not checked_any:
        print("FAIL: nothing was compared")
        return 1

    print("RANKS_BITWISE_IDENTICAL" if not failed else "RANKS_DIFFER")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
