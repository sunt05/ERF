# SLUCM single-column harness

Drives the SLUCM physics, `urban()`, offline for 48 h so that ERF's vendored
`module_sf_urban.F90` and the pristine upstream WRF `phys/module_sf_urban.F` it
was copied from can be run under identical inputs and their trajectories
compared.

This is item 6 of section 8 in
[`plan-slucm-coupling.md`](../../../Source/LandSurfaceModel/Noah-MP/dev/plan-slucm-coupling.md):
the regression guard on the vendoring and the kind handling.

## Result

At the shipped configuration the two builds are **bit-identical** — all 2880
steps, all twelve reported columns, 0 ULP:

```
steps compared : 2880
  TS     identical [gated]
  TR     identical [gated]
  ...
RESULT: BITWISE_IDENTICAL
```

The plan sets the gate at `1e-10` relative on `SH/LH/G/TS/TR/TB/TG`. That is the
*allowance*; the target is exact. Both builds compile the same Fortran with the
same promotion flags, so any non-zero difference would be a finding, and
`compare_trajectories.py` reports bitwise equality separately from the gate for
that reason.

## Why bit-identity is the right expectation here

The vendored copy carries 31 `! ERF:` edits against upstream. Four of them
change behaviour — they repair reads of undefined variables:

- `TGE` and `SROOT` are copied in `urban()`'s prologue before upstream ever
  assigns them.
- `tloc`/`tloc2` are computed only under `ahoption == 1` but read unconditionally
  by the `IRI_SCHEME == 1` irrigation block.
- `SW`/`BW` are read before assignment in a branch whose guarded and fallback
  bodies were identical.

ERF's `URBPARM.TBL` sets `AHOPTION 0`, `ALHOPTION 0`, `TREEOPTION 0`,
`GROPTION 0`, `IRI_SCHEME 0`, and `NoahmpUrbanInit` aborts if `TREEOPTION` or
`GROPTION` is 1. So in the configuration ERF actually runs, every one of those
edits touches code that is either not executed or whose result is not consumed —
hence identical bits. The comparison would *not* be well posed with the tree or
irrigation options on: there upstream reads undefined memory, and there is
nothing stable to compare against.

## Provenance of the reference

The vendored header originally recorded "no VCS metadata accompanied the copy".
The revision was recovered by md5-scanning all 49 revisions of
`phys/module_sf_urban.F` across every branch of `wrf-model/WRF`:

- `8fa379b458415a120fa8328cff78194b42c64a35` — 2026-05-19,
  *"Fixing a scheme-guard bug in urban NbS initialization (#2329)"*
- 5620 lines, md5 `ad16e13c8f24b30c9fa815fe10702db1`

`make fetch-upstream` retrieves exactly that blob and prints the md5 so the
reference can be re-established on any machine. The upstream file is never
edited; `wrf_stubs.F90` supplies the two modules it `use`s so it compiles
without a WRF tree.

## Running it

```bash
make fetch-upstream     # once; needs network
make run                # builds both sides, runs 48 h, compares
```

`make erf` builds the ERF side alone and works with no WRF source present.

Plain gfortran — no NetCDF, no MPI, no AMReX. ERF proper cannot be built on
macOS because `ERF_NCInterface.H` includes `netcdf_par.h` unconditionally and
needs a parallel netcdf-c; nothing here includes it, so a laptop is a fine host.
The case is single-node CPU work and should not be submitted to a GPU
allocation.

## Files

- `slucm_column.F90` — the driver. Mirrors `NoahmpUrbanDriverMainMod`'s call
  site so the state `urban()` sees here is the state it sees inside ERF, but
  calls `urban()` directly rather than through the Noah-MP IO layer, which would
  drag in the coupling this test holds constant.
- `forcing_48h.csv` — every `INTENT(IN)` field, 17 significant digits, so
  neither build can diverge through the driver. Regenerate with
  `make_forcing.py`.
- `golden_trajectory.csv` — committed baseline. The coming momentum and
  radiation fixes (section 4 of the plan) will change these numbers; the point
  of the baseline is that they show up as *intended* diffs rather than silent
  drift.
- `compare_trajectories.py` — bitwise comparison plus the gated relative check.

## Known open observation

Building either side with `-finit-real=snan -ffpe-trap=invalid,zero,overflow`
traps. **Both** sides do, ERF included, so this is not a demonstration of the
four upstream defects above; it indicates further reads of uninitialised locals
in `urban()` in this configuration. Whether they originate in the module or in
this driver has not been established, and no claim is made either way. The 3-D
smoke test cannot investigate this at all — Noah-MP writes through HDF5 and
`H5T__init_native_float_types` raises `FE_INVALID` during library init — so this
harness is the only place it can be chased.

## What this does not test

Physics. The trajectory is a synthetic mid-latitude summer day, not a site.
Separating "we broke the port" from "we broke the coupling" is item 5 of the
plan's section 8: the offline benchmark against the Grimmond et al. (2010, 2011)
urban energy-balance comparison forcing, judged against the published SLUCM
submission and the observations. `forcing_48h.csv`'s column layout is the one
that reader should produce.
