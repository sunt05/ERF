# SLUCM smoke test

A small end-to-end check that the WRF single-layer urban canopy model
(`SF_URBAN_PHYSICS = 1`) actually runs as an urban tile inside ERF's Noah-MP
land surface model, changes what it should change, and conserves energy while
doing it.

It is deliberately tiny — a 10 x 10 x 20 synthetic domain, 1800 steps of 1 s,
about 30 s of wall clock on one core — because its job is to catch a coupling
that has gone inert or gone wrong, not to validate urban physics. Validation
against observations is a separate exercise; see section 8 of
[`plan-slucm-coupling.md`](../../../Source/LandSurfaceModel/Noah-MP/dev/plan-slucm-coupling.md).

## Running it

Build ERF with Noah-MP, NetCDF and RRTMGP:

```bash
cmake -B build -DERF_ENABLE_MPI=ON -DERF_ENABLE_NETCDF=ON -DERF_ENABLE_NOAHMP=ON -DERF_ENABLE_RRTMGP=ON -DERF_ENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release
```

then

```bash
./run_smoke.sh --exe /path/to/build/Exec/erf_exec --work /somewhere/slucm_smoke
```

`run_smoke.sh` generates the input, stages the tables and the RRTMGP
coefficient files out of the submodules, runs the case twice, and calls
`check_smoke.py`. Nothing is downloaded and nothing outside the work directory
is written. `python3` with `netCDF4` is needed for the generator and the
checker; `--no-check` skips the latter.

On a machine whose MPI can only be started through the batch system, set
`ERF_MPI_LAUNCH` to a launcher that takes the rank count last:

```bash
ERF_MPI_LAUNCH="srun -n" ./run_smoke.sh --exe ... --work ...
```

This is needed on Cray with Slurm, where there is no usable `mpirun` and a bare
exec inside an allocation fails in PMI rather than running serially — so the
single-rank path needs redirecting too, not just `--np > 1`. Leaving the
variable unset reproduces the previous behaviour exactly.

The NetCDF stack must be a *parallel* build: `Source/IO/ERF_NCInterface.H`
includes `<netcdf_par.h>` unconditionally, and a serial netcdf-c does not
install that header. Ubuntu's `libnetcdf-dev` and Homebrew's `netcdf` are both
serial; `.github/workflows/gcc-noahmp.yml` builds the stack with Spack for
exactly this reason.

## What it runs

`Tests/Utils/make_synthetic_wrfinput.py` fabricates `wrfinput_urban_d01`: a
1 km grid, flat, USGS land use, a 4 x 4 block of urban cells
(`IVGTYP = ISURBAN`, `FRC_URB2D = 0.9`, `UTYPE = 2`) in an otherwise rural
domain, 1 m/s westerly, dated 2020-07-01 19:00 UTC — local solar noon at the
domain's longitude of -105, so the sun is high and the canyon geometry matters.
No `wrfbdy` is produced; with no boundary file ERF sets `use_real_bcs = false`,
which is what makes the periodic lateral boundaries in `inputs_slucm` legal.

The same `inputs_slucm` and the same `namelist.erf` drive both runs. The only
difference is `SF_URBAN_PHYSICS`, which `run_smoke.sh` rewrites to 0 for the
control, so the two configurations cannot drift apart in any other setting. The
control is *not* "no urban": it is Noah-MP's existing bulk-urban treatment,
which swaps in an urban roughness and urban soil parameters and stops there.
The test therefore measures what the canopy model adds over what ERF already
had.

## What it asserts

Read from the Noah-MP land output, `lnd<step>/Level_0.nc`:

1. **Both runs completed and every land field is finite** at every sample.
   This stands in for `amrex.fpe_trap_invalid`, which cannot be used here:
   Noah-MP writes through `NF90_NETCDF4 | NF90_MPIIO`, so HDF5 always
   initialises and `H5T__init_native_float_types` raises `FE_INVALID` while
   probing native float types, tripping the trap during library init.
2. **`FRC_URB2D` survives wrfinput -> ERF -> Noah-MP** unchanged, and the urban
   diagnostics appear in the urban run and only there.
3. **The SLUCM tile conserves energy**: `RN_URB = SH_URB + LH_URB - G_URB`
   (`urban()` defines `G = -FLXG*697.7*60`, hence the minus sign).
4. **Noah-MP's own budget still closes on the rural cells**,
   `RN = SH + LH + G + CANHS`, in both runs — so the urban work has not
   disturbed anything outside the urban tile.
5. **The urban cells respond**: surface temperature, ground flux, sensible
   heat, momentum stress and roughness all move away from the control.
6. **That response stays local**: rural roughness is bit-identical, and the
   rural flux response stays well below the urban one.

Note on check 3. The blend in `NoahmpUrbanDriverMainMod` overwrites
`HFX`/`LH`/`GRDFLX`/`TSK`/`EMISS`/`ALBEDO`/`TAU_*` with tile-weighted values but
leaves `FSAXY`/`FIRAXY`/`SAGXY`/`SAVXY` as the *rural tile's*. So
`FSAXY - FIRAXY` is not the blended cell's net radiation, and pairing it with
the blended fluxes on an urban cell leaves a residual of order
`FRC * (RN_rural - RN_urban)` — about 5% of RN in this case — that says nothing
about conservation. The check uses `RN_URB2D` and its companions instead, which
are `urban()`'s own budget terms.

Note on check 6. Rural cells are *not* expected to be bit-identical between the
two runs, and asserting that they are would be asserting wrong physics: the
urban block heats the air above it, that air advects, the pressure field couples
the domain within one acoustic substep, and RRTMGP sees a different surface
albedo. Only roughness is a pure per-column surface property, so only roughness
gets an exact-equality assertion.

## Running it on more than one rank

`--np N` works for any N that divides the 10-cell x extent. Two overrides go in
automatically, and both are about ERF and Noah-MP rather than about SLUCM:

- Noah-MP asserts that every box starts at `klo`
  (`ERF_NOAHMP_Init.cpp:241`, "z-decomposed grids are unsupported"), so the
  split has to be in x only, and `amr.blocking_factor` has to drop to 1 because
  10 cells do not divide by the default 8.
- **ERF's 2-D NetCDF plotfile writer hangs on more than one rank.** With
  `erf.plotfile2d_type_1 = netcdf` on 2 ranks this case stalls at the first
  `plt2d` write and never returns; the Noah-MP land NetCDF written just before
  it completes fine. Setting `erf.plot2d_int_1 = -1` lets the same run finish,
  so `run_smoke.sh` switches the 2-D plotfile off for `--np > 1`. It has not
  been isolated to a case without Noah-MP, so calling it a pre-existing ERF bug
  rather than an interaction is still an inference, not a finding.

With those in place, 1 rank and 2 ranks agree **bit for bit** on all 54 land
fields, which is the decomposition-independence check worth having. Run it with
`compare_ranks.py`, which takes two work directories and compares every field:

```bash
./run_smoke.sh --exe ... --work /tmp/np1 --np 1
./run_smoke.sh --exe ... --work /tmp/np2 --np 2
./compare_ranks.py /tmp/np1 /tmp/np2
```

Bit-for-bit is deliberately the bar: a decomposition must not perturb a
per-column land calculation at all, so a difference of any size is a bug in halo
exchange, in Noah-MP's box assumptions, or in the urban tile's indexing — not
something to absorb into a tolerance. NaN compares equal to NaN, since masked
land points are everywhere and would otherwise swamp the result.

## Reference output

On one core, gfortran 11.5, Linux, at step 1800:

```
SLUCM tile energy budget (RN_URB = SH_URB + LH_URB - G_URB)
  |RN_URB| max            714.40 W m^-2
  max residual              0.043 W m^-2

Noah-MP bulk budget on rural cells (RN = SH + LH + G + CANHS)
  urban run max residual    4.7e-05 W m^-2
  control  max residual     4.3e-05 W m^-2

urban response (max |urban - control|)       urban        rural     ratio
  TSK                                        18.57 K       1.71 K    10.8
  GRDFLX                                    290.9 W/m2    18.3       15.9
  HFX                                       407.4 W/m2    40.0       10.2
  TAU_EW                                      0.634 Pa     0.142      4.5
  ZNT                                         0.86 m       0.0        inf
```

The numbers will move with any change to the physics or the fixture; the
assertions are on orders of magnitude and on the two conservation identities,
not on these values.
