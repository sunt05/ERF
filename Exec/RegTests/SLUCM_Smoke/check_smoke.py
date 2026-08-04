#!/usr/bin/env python3
"""Check the SLUCM urban run against the bulk-urban control.

Reads the Noah-MP land NetCDF output (``lnd<step>/Level_0.nc``) written by
``NoahmpWriteLandMod.F90`` from both run directories laid down by
``run_smoke.sh`` and checks that

1. both runs produced the full set of land samples and every field is finite;
2. ``FRC_URB2D`` survives the whole path -- wrfinput -> ERF -> Noah-MP -- and
   the urban diagnostics appear only when ``SF_URBAN_PHYSICS = 1``;
3. the SLUCM tile closes its own surface energy budget,
   ``RN_URB = SH_URB + LH_URB - G_URB`` (``urban()`` defines
   ``G = -FLXG*697.7*60``, hence the minus sign);
4. Noah-MP's bulk budget still closes on the rural cells,
   ``RN = SH + LH + G + CANHS``, in both runs;
5. the urban cells respond to the urban tile -- surface temperature, ground
   flux, sensible heat, momentum stress and roughness all move away from the
   control;
6. that response stays localised on the urban cells -- rural roughness is
   bit-identical, and the rural flux response stays far below the urban one.

Exit status is 0 only if every assertion holds.

Two things are worth knowing before reading a residual out of this file.

*The blend does not cover the radiation diagnostics.*  ``NoahmpUrbanDriverMain``
overwrites ``HFX``/``LH``/``GRDFLX``/``TSK``/``EMISS``/``ALBEDO``/``TAU_*``
with tile-weighted values but leaves ``FSAXY``/``FIRAXY``/``SAGXY``/``SAVXY`` as
the *rural tile's*.  So ``FSAXY - FIRAXY`` is not the blended cell's net
radiation and combining it with the blended fluxes on an urban cell leaves a
residual of order ``FRC * (RN_rural - RN_urban)`` -- around 5% of RN in this
case -- which says nothing about whether SLUCM conserves energy.  That is why
check 3 uses ``RN_URB2D`` and friends, which are ``urban()``'s own budget terms.

*Index convention.*  The Fortran writer declares its variables ``(/nx, ny/)``
and the netCDF Fortran API reverses dimension order, so on disk -- and
therefore in numpy -- these arrays are ``[j, i]``.  The wrfinput fields are
written ``("Time", "south_north", "west_east")``, i.e. ``[0, j, i]``.  Both are
indexed the same way and no transpose is needed.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import numpy as np
from netCDF4 import Dataset

# Fields checked for finiteness in every sample of both runs.
FINITE_FIELDS = [
    "TSK", "HFX", "LH", "GRDFLX", "FIRAXY", "FSAXY", "SAGXY", "SAVXY",
    "EMISS", "TAU_EW", "TAU_NS", "Z0", "ZNT", "TSLB", "SMOIS", "TRADXY",
]

# Diagnostics that exist only under SF_URBAN_PHYSICS > 0.
URBAN_DIAGS = ["SH_URB2D", "LH_URB2D", "G_URB2D", "RN_URB2D", "TS_URB2D",
               "UST_URB2D", "FRC_URB2D", "UTYPE_URB2D"]

# (field, minimum |urban - control| required on urban cells, units)
# Thresholds sit far above round-off and far below "the physics blew up": the
# point is to catch a silently inert urban tile, not to pin numbers that a
# legitimate physics change would move.
URBAN_RESPONSE = [
    ("TSK",    1.0e-3, "K"),
    ("GRDFLX", 1.0e-1, "W m^-2"),
    ("HFX",    1.0e-1, "W m^-2"),
    ("TAU_EW", 1.0e-4, "Pa"),
    ("ZNT",    1.0e-4, "m"),
]

# Fields the locality check is run on.  LH is deliberately absent: the SLUCM
# tile returns LH_URB = 0 on a dry canyon and the bulk-urban control returns
# LH = 0 outright, so the urban-cell change in LH is only the 1 - FRC rural
# remnant (~20 W m^-2 here) -- the same order as the rural cells' own
# atmospheric response.  It carries no locality signal in a dry case, and
# asserting on it would just be asserting on the moisture state of the fixture.
LOCALITY_FIELDS = ["TSK", "HFX", "GRDFLX", "TAU_EW", "TAU_NS", "ZNT"]

# Absolute tolerances, W m^-2.  Both budgets are evaluated in single precision
# (the land file is NF90_FLOAT) on terms of order 1000, so ~1e-1 is the
# representation floor, not a slack allowance.
TILE_CLOSURE_TOL = 1.0
BULK_CLOSURE_TOL = 1.0

# How much larger the urban-cell response must be than the rural-cell response.
# The rural response is real physics (see check 6), so this is a contrast
# requirement, not a tolerance.
LOCALITY_RATIO = 3.0


def land_samples(rundir: str) -> list[tuple[int, str]]:
    """Return [(step, path)] for every lnd<step>/Level_0.nc, ordered by step."""
    out = []
    for d in glob.glob(os.path.join(rundir, "lnd*")):
        m = re.match(r"lnd(\d+)$", os.path.basename(d))
        f = os.path.join(d, "Level_0.nc")
        if m and os.path.isfile(f):
            out.append((int(m.group(1)), f))
    return sorted(out)


def read(path: str, names=None) -> dict:
    with Dataset(path) as nc:
        keys = list(nc.variables) if names is None else \
            [n for n in names if n in nc.variables]
        return {n: np.array(nc.variables[n][:], dtype=np.float64) for n in keys}


def variables(path: str) -> list[str]:
    with Dataset(path) as nc:
        return list(nc.variables)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work", required=True,
                   help="directory laid down by run_smoke.sh")
    args = p.parse_args(argv)

    urban_dir = os.path.join(args.work, "urban")
    rural_dir = os.path.join(args.work, "rural")
    wrfinput = os.path.join(args.work, "wrfinput_urban_d01")

    failures: list[str] = []

    with Dataset(wrfinput) as nc:
        if "FRC_URB2D" not in nc.variables:
            print("FAIL: the wrfinput carries no FRC_URB2D")
            return 1
        frc_in = np.array(nc.variables["FRC_URB2D"][0, :, :], dtype=np.float64)
    mask = frc_in > 0.0
    n_urb, n_rur = int(mask.sum()), int((~mask).sum())
    print(f"urban cells: {n_urb} (FRC_URB2D = {frc_in[mask].max():.3f}), rural cells: {n_rur}")
    if n_urb == 0 or n_rur == 0:
        print("FAIL: the domain must contain both urban and rural cells")
        return 1

    su, sr = land_samples(urban_dir), land_samples(rural_dir)
    print(f"land samples: urban {[s for s, _ in su]}, control {[s for s, _ in sr]}")
    if not su or not sr:
        print("FAIL: one of the runs produced no land output at all")
        return 1
    if [s for s, _ in su] != [s for s, _ in sr]:
        failures.append("the two runs wrote land output at different steps")
    if len(su) < 3:
        failures.append(f"only {len(su)} land sample(s); the run stopped early")

    # ---- 1. finiteness -----------------------------------------------------
    for tag, samples in (("urban", su), ("control", sr)):
        for step, path in samples:
            for name, arr in read(path, FINITE_FIELDS).items():
                if not np.all(np.isfinite(arr)):
                    bad = int((~np.isfinite(arr)).sum())
                    failures.append(f"{tag} step {step}: {name} has {bad} non-finite values")

    last_u, last_r = su[-1][1], sr[-1][1]

    # ---- 2. the urban tile is switched on, and only in the urban run -------
    have_u = variables(last_u)
    missing = [v for v in URBAN_DIAGS if v not in have_u]
    if missing:
        failures.append(f"urban run land output lacks {', '.join(missing)}; "
                        f"either SF_URBAN_PHYSICS never reached Noah-MP or the "
                        f"writer's guard is wrong")
    leaked = [v for v in URBAN_DIAGS if v in variables(last_r)]
    if leaked:
        failures.append(f"control run wrote urban diagnostics ({', '.join(leaked)}) "
                        f"despite SF_URBAN_PHYSICS = 0")

    if not missing:
        du = read(last_u)
        frc_out = du["FRC_URB2D"]
        dfrc = float(np.abs(frc_out - frc_in).max())
        print(f"FRC_URB2D wrfinput -> Noah-MP: max |difference| {dfrc:.3g}")
        if dfrc > 1.0e-6:
            failures.append(f"FRC_URB2D reaching Noah-MP differs from the wrfinput "
                            f"by up to {dfrc:g}")
        utypes = np.unique(du["UTYPE_URB2D"][mask])
        print(f"urban type on urban cells: {utypes.astype(int).tolist()}")

        # ---- 3. the SLUCM tile closes its own budget -----------------------
        rn, sh, lh, g = (du["RN_URB2D"], du["SH_URB2D"],
                         du["LH_URB2D"], du["G_URB2D"])
        resid = np.abs(rn - (sh + lh - g))[mask]
        rmax = float(resid.max())
        print(f"\nSLUCM tile energy budget at step {su[-1][0]} "
              f"(RN_URB = SH_URB + LH_URB - G_URB):")
        print(f"  |RN_URB| max                  {float(np.abs(rn[mask]).max()):12.4f} W m^-2")
        print(f"  SH_URB  range                 [{float(sh[mask].min()):.2f}, {float(sh[mask].max()):.2f}] W m^-2")
        print(f"  LH_URB  range                 [{float(lh[mask].min()):.2f}, {float(lh[mask].max()):.2f}] W m^-2")
        print(f"  G_URB   range                 [{float(g[mask].min()):.2f}, {float(g[mask].max()):.2f}] W m^-2")
        print(f"  max residual                  {rmax:12.6f} W m^-2   "
              f"{'ok' if rmax <= TILE_CLOSURE_TOL else 'FAIL'}")
        if rmax > TILE_CLOSURE_TOL:
            failures.append(f"the SLUCM tile does not conserve energy: max residual "
                            f"{rmax:g} W m^-2 > {TILE_CLOSURE_TOL:g}")

    # ---- 4. Noah-MP's bulk budget still closes on the rural cells ----------
    print("\nNoah-MP bulk budget on rural cells (RN = SH + LH + G + CANHS):")
    for tag, path in (("urban run", last_u), ("control ", last_r)):
        e = read(path, ["FSAXY", "FIRAXY", "HFX", "LH", "GRDFLX", "CANHSXY"])
        if len(e) < 6:
            failures.append(f"{tag}: land output lacks the bulk-budget fields")
            continue
        resid = (e["FSAXY"] - e["FIRAXY"]) - (e["HFX"] + e["LH"] + e["GRDFLX"] + e["CANHSXY"])
        rmax = float(np.abs(resid[~mask]).max())
        print(f"  {tag} max residual         {rmax:12.6f} W m^-2   "
              f"{'ok' if rmax <= BULK_CLOSURE_TOL else 'FAIL'}")
        if rmax > BULK_CLOSURE_TOL:
            failures.append(f"{tag}: Noah-MP's own budget does not close on the rural "
                            f"cells (max residual {rmax:g} W m^-2); the urban work has "
                            f"broken something outside the urban tile")

    # ---- 5. urban cells respond -------------------------------------------
    print(f"\nurban response at step {su[-1][0]} (max |urban - control|):")
    du = read(last_u, [f for f, _, _ in URBAN_RESPONSE])
    dr = read(last_r, [f for f, _, _ in URBAN_RESPONSE])
    for name, thresh, unit in URBAN_RESPONSE:
        if name not in du or name not in dr:
            failures.append(f"{name} absent from the land output")
            continue
        diff = np.abs(du[name] - dr[name])
        mx_u, mx_r = float(diff[mask].max()), float(diff[~mask].max())
        ok = mx_u > thresh
        print(f"  {name:8s} urban {mx_u:12.6g} {unit:9s} rural {mx_r:12.6g}   "
              f"{'ok' if ok else f'FAIL (<= {thresh:g})'}")
        if not ok:
            failures.append(f"{name} does not respond to SF_URBAN_PHYSICS = 1 "
                            f"(max |diff| on urban cells {mx_u:g} <= {thresh:g})")

    # ---- 6. the urban signal stays where the urban cells are ---------------
    # Rural cells are *not* expected to be bit-identical between the two runs.
    # The urban block heats the air above it, that air advects, the pressure
    # field couples the whole domain within one acoustic substep, and RRTMGP
    # sees a different surface albedo -- so the rural surface fluxes respond
    # too.  Requiring identity there would be requiring wrong physics.
    #
    # What can be asserted is locality of two different strengths:
    #
    #   * ZNT is a pure per-column surface property with no atmospheric input,
    #     so a rural cell's roughness must be bit-identical.  Any change there
    #     is the urban tile writing outside its own cells.
    #   * the flux and temperature responses must remain concentrated on the
    #     urban cells, by a wide margin.  A blend that had leaked its weights
    #     would flatten that contrast.
    print(f"\nlocality at step {su[-1][0]} (max |urban - control|, urban vs rural cells):")
    lu = read(last_u, LOCALITY_FIELDS)
    lr = read(last_r, LOCALITY_FIELDS)
    for name in LOCALITY_FIELDS:
        if name not in lu or name not in lr:
            continue
        diff = np.abs(lu[name] - lr[name])
        mx_u, mx_r = float(diff[mask].max()), float(diff[~mask].max())
        if name == "ZNT":
            ok = mx_r == 0.0
            note = "must be exactly 0 off the urban cells"
        else:
            ok = mx_r <= mx_u / LOCALITY_RATIO
            note = f"rural must be <= urban / {LOCALITY_RATIO:g}"
        ratio = (mx_u / mx_r) if mx_r > 0 else float("inf")
        print(f"  {name:8s} urban {mx_u:12.6g}  rural {mx_r:12.6g}  ratio {ratio:9.3g}   "
              f"{'ok' if ok else 'FAIL'}")
        if not ok:
            failures.append(f"{name}: the urban signal is not localised ({note}; "
                            f"urban {mx_u:g}, rural {mx_r:g})")

    print()
    if failures:
        for f in failures:
            print(f"FAIL: {f}")
        return 1
    print("PASS: SLUCM runs inside ERF, conserves energy on the urban tile, keeps "
          "its response localised there, and leaves Noah-MP's rural budget intact.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
