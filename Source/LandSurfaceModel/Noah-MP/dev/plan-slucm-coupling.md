# Plan: couple the Single-Layer Urban Canopy Model (SLUCM) into ERF

> Status: proposal, not yet implemented · Owns: the design and phasing for adding
> a WRF-style single-layer urban canopy model as an **urban tile** alongside the
> Noah-MP rural tile · Companions: [`spec-noahmp-api.md`](spec-noahmp-api.md)
> (the `NOAHMP` class and coupling fields), [`spec-noahmp-gpu.md`](spec-noahmp-gpu.md)
> (per-step staging), `Submodules/Noah-MP/drivers/erf/dev/spec-add-coupled-variable.md`
> (how a variable crosses the C++ ↔ Fortran boundary).
> Public docs: <https://erf.readthedocs.io/en/latest/CouplingToNoahMP.html>

## 1. Context — what we have and what is missing

ERF's real-data land surface today is Noah-MP, selected with
`erf.land_surface_model = "NOAHMP"` (`Source/ERF.cpp:2469`) and driven by
`Source/LandSurfaceModel/Noah-MP/`. Over an urban grid cell it runs in **bulk-urban**
mode: `ConfigVarInTransferMod.F90:155-165` sets `VegType = ISURBAN_TABLE` and
`FlagUrban = .true.`, which only swaps in an urban roughness length
(`src/GroundRoughnessPropertyMod.F90:70-72`) and urban soil/thermal parameters.
There is no canyon geometry, no separate roof/wall/road energy balance, no urban
thermal storage, and no anthropogenic heat.

We want the WRF single-layer urban canopy model (Kusaka et al. 2001; Kusaka and
Kimura 2004; Chen et al. 2011 — WRF `sf_urban_physics = 1`) run as a **tile**
over the urban fraction of each cell, blended with the Noah-MP rural tile by
`FRC_URB2D`. That is the "mosaic" formulation:

```
HFX = FRC_URB2D · SH_URB + (1 − FRC_URB2D) · HFX_rural            [W m-2]
QFX = FRC_URB2D · LH_KINEMATIC_URB + (1 − FRC_URB2D) · QFX_rural  [kg m-2 s-1]
TSK / ALBEDO / EMISS / GRDFLX / QSFC / UST  blended the same way
```

Noah-MP supplies the rural side of that blend by being run with the *natural*
vegetation type over the whole cell, not the urban one.

## 2. Upstream survey — are we first?

Yes, for a mesoscale urban **parameterization**. Findings from `erf-model/ERF`
issues and pull requests (searched `urban`, `SLUCM`, `UCM`, `LCZ`, `canopy`,
`building`, `land surface`):

- **No urban canopy parameterization exists and none is proposed.** No open or
  closed issue/PR mentions SLUCM, BEP, BEP+BEM, TEB, LCZ, or `FRC_URB2D`.
- **ERF's urban story upstream is building-resolving, not parameterized.**
  Immersed forcing (`erf.buildings_type = "ImmersedForcing"`) plus embedded
  boundaries is actively developed: #2763 (IF terrain and buildings pathways),
  #2953 (immersed-boundary temperature forcing), #3464/#3468 (zero surface
  fluxes inside immersed cells), #3478 (wall model for buildings), #3482
  (implicit IF formulation, open draft). `Docs/sphinx_doc/ERFvsWRF.rst:69` is the
  only place the word "urban" appears in the docs, and it points at exactly this:
  *"Terrain and urban geometries may be simulated with immersed forcing or
  embedded (immersed) boundary techniques."*
- **The nearest precedent for a canopy parameterization is the forest drag model**
  (#1980, ported from AMR-Wind): `Source/SourceTerms/ERF_ForestDrag.{H,cpp}`
  builds an `a·Cd` field applied as a momentum sink in
  `ERF_MakeMomSources.cpp:747-783`. It has no energy budget — it is drag only.
- **All 2026 LSM effort upstream is Noah-MP hardening**, not new schemes: full
  two-way coupling (#3414/#3412), checkpoint/restart (#3255, #3317), sea-ice
  guard (#3286), coarse→fine flux interpolation (#3334), stress unit fix (#3428),
  expanded 2-D land diagnostics (#3465).

Two consequences for this plan:

1. There is no competing branch to coordinate with — but also no reviewer with
   context, so the design has to be self-justifying and testable against WRF.
2. ERF already "owns" the microscale urban problem. **SLUCM is the mesoscale and
   grey-zone complement, not a replacement**, and the plan must forbid running
   both over the same cells (§7).

On the WRF side, note that the tile blending no longer lives in WRF itself: WRF
master has moved Noah-MP into the `NCAR/noahmp` submodule, and the urban tile
driver now lives *inside that submodule* as
`drivers/wrf/NoahmpUrbanDriverMainMod.F90`. Only the SLUCM physics
(`phys/module_sf_urban.F`) and `run/URBPARM.TBL` remain in WRF proper. That
matters, because our Noah-MP submodule is a fork of the same repository.

## 3. What already exists in our tree

`Submodules/Noah-MP` is pinned to `erf-model/noahmp` @ `e0aed20`, a fork of
`NCAR/noahmp`. **Most of the Fortran-side urban plumbing is already there**,
inherited from the NCAR driver and never exercised because ERF passes
`SF_URBAN_PHYSICS = 0`:

| Already present in `drivers/erf/` | What it gives us |
|---|---|
| `NoahmpIOVarType.F90-mc` (133 urban references) | The complete SLUCM array set — `tr/tb/tg/tc/qc/uc_urb2d`, `xxx{r,b,g,c}_urb2d`, `trl/tbl/tgl_urb3d`, `sh/lh/g/rn/ts_urb2d`, `psim/psih/u10/v10/GZ1OZ0/AKMS/th2/q2/ust_urb2d`, `cmcr/tgr/tgrl/smr/drel*/flxhum*`, `frc_urb2d`, `utype_urb2d`, `mh/stdh/lp/lb/hgt/lf_urb2d` — plus the BEP/BEM and green-roof/PV sets, `sf_urban_physics`, `num_urban_*`, the derived `urban_map_*`, `ISURBAN`, `URBTYPE_beg`, `LCZ_1..11_TABLE`, `IRI_URBAN`, `GMT`, `JULDAY`, `HRANG`, `DECLIN` |
| `NoahmpIOVarInitMod.F90-mc` (187 urban references) | Allocation of all of the above, already gated on `SF_URBAN_PHYSICS == 1` vs `2/3` |
| `NoahmpReadNamelistMod.F90` | Reads `sf_urban_physics`, `use_wudapt_lcz`, `num_urban_hi`, `urban_atmosphere_thickness`; derives `urban_map_*` |
| `NoahmpReadTableMod.F90` + `NoahmpTable.TBL` | `ISURBAN_TABLE`, `NATURAL_TABLE`, `LCZ_1..11_TABLE` (`URBTYPE_beg = 50`, LCZ 51–61) |
| **`ConfigVarInTransferMod.F90:155-165`** | **The mosaic contract is already implemented.** For an urban cell with `SF_URBAN_PHYSICS > 0` it sets `VegType = NATURAL_TABLE` and `GVFMAX = 96 %` — i.e. Noah-MP runs the *rural tile*, and the comment says outright that "urban is handled by explicit urban scheme outside Noah-MP". `SF_URBAN_PHYSICS == 0` is today's bulk-urban path. |
| `NoahmpInitMainMod.F90:175-207`, `NoahmpDriverMainMod.F90:200-215` | Urban-point branching for LAI/biomass initialization and urban irrigation |
| `NoahmpReadLandMod.F90` | Reads the wrfinput NetCDF directly (`XLAT`, `XLONG`, `IVGTYP`, `VEGFRA`, `SHDMIN/MAX`, …) and the global attribute `ISURBAN` — the natural place to add urban static fields |
| `Exec/RegTests/WPS_Test/namelist.erf:56-57` | Already carries `SF_URBAN_PHYSICS = 0` and `USE_WUDAPT_LCZ = 0` |

**Exactly three Fortran pieces are missing:**

1. **`module_sf_urban.F`** — the SLUCM physics. Not in the noahmp repo; it is WRF
   `phys/module_sf_urban.F` (5620 lines). Its only external dependencies are
   `wrf_error_fatal`/`wrf_message` (via a preprocessor shim it already defines at
   the top of the file) and `piconst` from `module_model_constants` — both
   trivially satisfied. It also brings `urban_param_init` (reads `URBPARM.TBL`)
   and `urban_var_init` (cold-start of the urban state), which WRF calls from
   `module_physics_init.F` and which we must call from `NoahmpInitMain`.
2. **`NoahmpUrbanDriverMainMod.F90`** — present in the fork at `drivers/wrf/`
   (839 lines), absent from `drivers/erf/`. Its `sf_urban_physics == 1` branch
   (lines 411–655) *is* the mosaic blender, line-for-line reusable.
3. **`URBPARM.TBL`** — shipped next to `NoahmpTable.TBL` in the run directory.

## 4. Gaps on the ERF (C++) side

| Gap | Detail |
|---|---|
| **Coupled variables** | The urban *state* arrays never need to reach C++ — they are internal Noah-MP state. What must cross: **in** — `FRC_URB2D`, `UTYPE_URB2D`, the NUDAPT morphology (`LP/LB/HGT/MH/STDH/LF_URB2D`), `XLONG`, and the time scalars `JULIAN`, `YR`, `GMT`, `DECLIN`; **out** — nothing new is strictly required, because SLUCM blends in place into `HFX`, `LH`, `TAU_EW/NS`, `TSK`, `EMISS`, `GRDFLX`, `ALBEDO`, all already coupled. Optional diagnostics: `TS_URB2D`, `SH/LH/G/RN_URB2D`, `TR/TB/TG/TC_URB2D`, `UST_URB2D`. Each is one line in the `@NoahmpMacro:Source m_noahmpio { … }` block of `NoahmpIO.H-mc`. |
| **Static urban input** | `FRC_URB2D`, `URB_PARAM`, `LU_INDEX`, and LCZ are read **nowhere** in ERF (`ERF_ReadFromWRFInput.cpp:125-129` reads `IVGTYP`/`ISLTYP` only; metgrid reads no land-use at all). |
| **`urb_frac_lev`** | Already declared (`Source/ERF.H:1033`), allocated and `setVal(1)` (`Source/ERF_MakeNewArrays.cpp:501-504`), resized (`ERF_Constructors.cpp:300`) — and **never filled or read anywhere**. It is a ready-made slot for `FRC_URB2D`. |
| **Solar geometry** | SLUCM needs `DECLIN_URB`, `COSZ_URB2D`, `OMG_URB2D` (hour angle) and `XLAT_URB2D`. ERF already couples `cos_zenith_angle` (`LsmData_NOAHMP`) and has `orbital_decl` / `orbital_cos_zenith` in `Source/PhysicsInterfaces/Radiation/ERF_OrbCosZenith.H` — declination and hour angle come straight out of that machinery. |
| **Roughness / displacement height** | The real structural gap — but smaller than it looks. `Z0` and `ZNT` are **already coupled members** of the `@NoahmpMacro:Source` block (`NoahmpIO.H-mc:169-170`), and Noah-MP already fills `ZNT`; SLUCM blends into it in place. What is missing is entirely on ERF's side: `ERF_NOAHMP_Fields.H` never reads `ZNT` into an `LsmData` field, `SurfaceLayer` takes `z0` only from `erf.most.z0` or `erf.most.roughness_file_name`, and `rough_type_land` is **hard-`Abort`ed** unless `"constant"` (`ERF_SurfaceLayer.H:205-215`). Separately, there is **no displacement height `d0` anywhere in ERF** — the log law is always `log(zref/z0)` (`ERF_MOSTStress.H:174`). |
| **Advance dispatch** | `Source/TimeIntegration/ERF_AdvanceLSM.cpp:13` hard-codes `LandSurfaceType::NOAHMP` for the `Advance_With_State` path. Fine if SLUCM lives inside the Noah-MP model; a separate `LandSurfaceType` would have to be added here. |

> **The roughness gap is not a Phase-1 blocker.** SLUCM's effect reaches the
> atmosphere entirely through `t_flux`, `q_flux`, `tau13`, `tau23`, and
> `SurfaceLayer::compute_sfc_params_from_lsm_fluxes` (`ERF_SurfaceLayer.cpp:1016-1082`)
> back-derives `u*`, `θ*`, `q*`, and `L` from those. So momentum and heat coupling
> are already correct without touching `z0`. `z0` matters only for the PBL schemes'
> `get_z0` and for the MOST fallback in cells the LSM did not process.

## 5. Architecture decision

**Recommendation: SLUCM as an urban tile inside the Noah-MP Fortran driver**
(design A), mirroring `NCAR/noahmp`'s `drivers/wrf/NoahmpUrbanDriverMainMod.F90`.

The alternative (design B) is an ERF-native C++ scheme — a new
`LandSurfaceType::SLUCM` deriving from `NullSurf`. Reject it, because:

- The mosaic contract already exists in `ConfigVarInTransferMod.F90:155-165`.
  Design B would have to re-implement "Noah-MP over the rural fraction" from
  outside the Fortran model.
- `LandSurface` holds **one model per level**, and `Set_Lev0_*_Ptr` aborts if the
  model type differs across levels (`ERF_LandSurface.H:153-185`). You cannot run
  "Noah-MP for the rural tile plus SLUCM for the urban tile" as two `NullSurf`
  instances.
- `TSK`, `ALBEDO`, and `EMISS` feed RRTMGP. Blending inside the Fortran driver
  keeps one source of truth; design B needs a second blending layer in ERF.
- Bit-level traceability to WRF for validation — the whole point of adopting a
  published scheme rather than writing a new one.

Costs of A, accepted knowingly: we vendor WRF Fortran into the `erf-model/noahmp`
fork (WRF is in the public domain, so this is clean, but provenance and the
upstream revision must be recorded in the file header), and SLUCM stays host-side
— which matches how Noah-MP already runs, since `Advance_With_State` already does
a device→host staging round-trip (`spec-noahmp-gpu.md`). The `@internal` tier
described in `spec-add-coupled-variable.md` §"Step 0" is designed for exactly this
case: urban state arrays get generated storage, allocation, and (later) device
residency at **zero ABI cost**.

**Where design B *would* be right**: a multi-layer scheme (BEP / BEP+BEM) that
injects drag, heat, and TKE at several model levels is not a pure LSM tile — it
reaches into `ERF_MakeMomSources.cpp` the way `ForestDrag` does. That is future
work; the `a_u_bep`, `a_t_bep`, `sf_bep`, `vl_bep` arrays are already declared in
`NoahmpIOVarType.F90-mc`, so the door is open.

## 6. Phased plan

### Phase 0 — groundwork (no code)

- Record the WRF revision SLUCM is taken from, and open a tracking issue on
  `erf-model/ERF` describing the mosaic design so upstream sees it before the PR.
- Decide the exchange-coefficient question (§7) and the reference-height guard
  (§7) — both change the Fortran argument list, so settle them first.
- Build a WRF `sf_urban_physics = 1` reference case over the same domain as
  `Exec/RegTests/WPS_Test` to validate against.

### Phase 1 — Fortran: vendor SLUCM into the submodule fork

Work in `erf-model/noahmp`, branch off the pinned `e0aed20`:

1. Add `drivers/erf/module_sf_urban.F` from WRF `phys/`, with:
   - a small `NoahmpUrbanShimMod.F90` (or reuse `NoahmpFatalMod.F90`) providing
     `wrf_message`/`wrf_error_fatal` → `NoahmpIO_abort`, and `piconst`;
   - **explicit kind conversion**: the submodule builds with `DOUBLE_PREC` and
     `c_kind_noahmp`, while `module_sf_urban.F` uses bare `REAL`. Convert
     declarations to `real(kind=kind_noahmp)` rather than relying on
     `-fdefault-real-8`, which is not portable across the compilers ERF supports.
2. Add `drivers/erf/NoahmpUrbanDriverMainMod.F90`, copied from `drivers/wrf/` and
   reduced to the `sf_urban_physics == 1` branch (lines 411–655 plus the
   `*_RURAL` snapshot at 400–409 and the combined-radiation epilogue). Keep the
   BEP/BEM branches out of Phase 1; leave a clearly marked stub.
3. Call it from `drivers/erf/NoahmpDriverMainMod.F90` **after** the `JLOOP`/`ILOOP`
   over land columns, guarded by `if (NoahmpIO%SF_URBAN_PHYSICS == 1)`.
4. Call `urban_param_init` and `urban_var_init` from `NoahmpInitMainMod.F90`
   under the same guard; ship `URBPARM.TBL` next to `NoahmpTable.TBL`.
5. Extend `NoahmpReadLandMod.F90` to read `FRC_URB2D`, `URB_PARAM` (or the
   individual `LP/LB/HGT/MH/STDH/LF_URB2D` fields), and `LU_INDEX`/LCZ from the
   wrfinput file, `NOT_FATAL` with sane defaults so non-urban runs are unaffected.
6. Add the new files to `drivers/erf/CMakeLists.txt` and `drivers/erf/Makefile`.

### Phase 2 — interface: coupled variables

Follow `Submodules/Noah-MP/drivers/erf/dev/spec-add-coupled-variable.md` exactly.

- Add one line per crossing variable to the `@NoahmpMacro:Source m_noahmpio { … }`
  block in `NoahmpIO.H-mc` — inputs first (`FRC_URB2D`, `UTYPE_URB2D`, morphology,
  `XLONG`), then the scalars (`JULIAN`, `YR`, `GMT`, `DECLIN`), then any
  diagnostics we choose to export.
- Because `FRC_URB2D` and friends **already exist as hand-written Fortran
  allocatables**, this is the *promotion* path: delete the hand-written
  declaration in `NoahmpIOVarType.F90-mc` and the hand-written `allocate()` in
  `NoahmpIOVarInitMod.F90-mc` for each promoted variable, or the build fails with
  a duplicate-component error.
- `make codegen && make codegen-check` must be clean.
- Sentinel-guard the namelist scalars (`if (NoahmpIO%x == undefined_real) …`) so a
  C++-supplied value wins over the namelist.
- Bump the `Submodules/Noah-MP` pin in ERF.

### Phase 3 — ERF driver

`Source/LandSurfaceModel/Noah-MP/`:

- `ERF_NOAHMP_Fields.H` — add the new forcing rows to `NOAHMP_INPUT_2D_FIELDS`
  and, if we export urban diagnostics, rows to `NOAHMP_LSMDATA_FIELDS`,
  `NOAHMP_OUTPUT_2D_FIELDS_TAIL`, and `NOAHMP_RESULT_FIELDS`. **Append only** —
  enum order is an invariant.
- `ERF_NOAHMP_Advance.cpp` — in `stage_forcing`, compute the solar geometry
  (`orbital_decl` for `DECLIN`, the hour angle for `OMG`) and stage `JULIAN`,
  `YR`, `GMT`. These are per-step scalars, so set them next to `itimestep` rather
  than through the pinned FAB.
- `ERF_NOAHMP_Init.cpp` — after `ReadLandMain()`, copy `FRC_URB2D` into
  `urb_frac_lev[lev][0]` so the rest of ERF can see the urban fraction, and
  assert that `SF_URBAN_PHYSICS == 1` implies a per-level land file (see the AMR
  risk in §7).
- `Source/ERF_MakeNewArrays.cpp` / `Source/ERF.H` — no new fields needed;
  `urb_frac_lev` just stops being dead.
- Optionally register the urban diagnostics in the 2-D plotfile catalog
  (`Source/IO/ERF_Plotfile2DCatalog.cpp`), following the pattern added by
  upstream #3465.

### Phase 4 — roughness and displacement height (optional, structural)

Only worth doing once Phase 1–3 validate. Two independently sized pieces:

- **Roughness (cheap).** `ZNT` already crosses the boundary, so this is
  ERF-side only: add a `znt` row to `NOAHMP_LSMDATA_FIELDS` /
  `NOAHMP_OUTPUT_2D_FIELDS_TAIL` / `NOAHMP_RESULT_FIELDS`, add
  `RoughCalcType::LSM` and accept `erf.most.roughness_type_land = "lsm"`
  (`ERF_SurfaceLayer.H:205-215`), then fill `z_0[lev]` from that field the way
  `get_lsm_tsurf` fills `t_surf`. This is useful on its own even without SLUCM,
  since Noah-MP already varies `ZNT` by land-use type.
- **Displacement height (expensive).** Introduce a `d0` field and use
  `log((zref − d0)/z0)` in `ERF_MOSTStress.H`. This touches every flux iterator,
  so it needs its own spec and its own PR — **do not fold it into the SLUCM PR**.

### Phase 5 — tests, docs, regression case

- **Unit test**, mirroring `Tests/Unit/LandSurfaceModel/Noah-MP/ERF_GTestNoahMPResultPolicy.cpp`:
  extract the mosaic blend into a small pure header (e.g.
  `ERF_NOAHMP_UrbanTile.H`, a `blend(frc, urban, rural)` helper) and test it
  directly — `frc = 0` reproduces the rural value bitwise, `frc = 1` the urban
  value, sentinel propagation is correct. Register it in
  `Tests/Unit/CMakeLists.txt` inside the existing `if(ERF_ENABLE_NOAHMP)` block.
- **Fortran unit test** under `Submodules/Noah-MP/drivers/erf/tests/` for
  `urban_param_init` reading `URBPARM.TBL`.
- **Regression case**: a `WPS_Test`-derived urban case with
  `SF_URBAN_PHYSICS = 1`, added to `Tests/CTestList.cmake`. Note there is
  currently **no** Noah-MP regression test and no CI job setting
  `ERF_ENABLE_NOAHMP` — adding one is a prerequisite, and is worth doing on its
  own merits.
- **Docs**: `Docs/sphinx_doc/CouplingToNoahMP.rst` gains an "Urban canopy" section
  (build flags, `namelist.erf` keys, required wrfinput fields, `URBPARM.TBL`);
  `Docs/sphinx_doc/Inputs.rst:1918-1938` gains the urban rows;
  `Docs/sphinx_doc/ERFvsWRF.rst:69` is updated so urban is no longer described as
  immersed-forcing-only. Add `spec-slucm-tile.md` next to this plan and index it
  in `dev/README.md`.

## 7. Risks and open questions

1. **Double counting with immersed forcing.** If `erf.buildings_type = "ImmersedForcing"`
   and SLUCM are both active on the same level, buildings are represented twice —
   once as resolved geometry (`lmask == 2`, zero surface flux, `ERF_SurfaceLayer.cpp:591,632,692,735`)
   and once as a canopy parameterization. **Abort at parameter-check time.**
2. **Reference height in the grey zone.** SLUCM assumes the forcing level is above
   the canopy. The ERF driver currently sets `DZ8W = 2·ZLVL` and
   `P8W(:,2,:) = P8W(:,1,:)` (`NoahmpDriverMainMod.F90:52-66`), and SLUCM takes
   `ZA_URB = 0.5·DZ8W(i,1,j)`. With typical ERF vertical resolution the first
   level can sit *inside* the canopy. Add an explicit guard comparing `ZA_URB`
   against `HGT_URB2D` and either abort or document a minimum first-level height.
3. **Exchange coefficients.** WRF passes `CHS`, `CHS2`, `CQS2` from the surface-layer
   scheme into SLUCM; ERF's Noah-MP driver does not plumb them. Either derive them
   from `SurfaceLayer`'s `u*`/`t*` or let SLUCM use its internal `SFCDIF_URB` /
   `mos` / `louis` options (`AHOPTION`/`SFCDIF` in `URBPARM.TBL`). **Decide in
   Phase 0** — it changes the argument list.
4. **Precision.** See Phase 1.1. Getting this wrong produces silently wrong
   numbers, not a crash.
5. **AMR.** `interp_from_lev0` (`ERF_NOAHMP_Advance.cpp:23-52`) interpolates all
   LSM data and fluxes from level 0 when a fine level has no NetCDF land file.
   Urban morphology is strongly heterogeneous and interpolating `FRC_URB2D` across
   levels is not defensible. **Require a per-level land file when SLUCM is on.**
6. **GPU.** SLUCM is host-only. No new sync points are introduced (the staging
   round-trip already exists), but the GPU-offload plan
   (`Submodules/Noah-MP/drivers/erf/dev/plan-cpp-interface.md`) must account for it.
7. **Vendoring drift.** Once `module_sf_urban.F` lives in our fork it will drift
   from WRF. Record the source revision in the file header and keep the file
   otherwise unmodified apart from the kind and shim changes, so a future
   re-sync is a readable diff.
8. **Licensing.** `erf-model/noahmp` carries the **UCAR Noah-MP license** (royalty-free,
   no resale, attribution required on derived works). WRF is public domain, so
   vendoring `module_sf_urban.F` and `URBPARM.TBL` into the fork is clean — but
   add the WRF attribution header and cite Kusaka et al. (2001) / Chen et al.
   (2011) in the docs, as the UCAR license requires for derived works.
9. **Build prerequisites.** The submodule's CMake runs `tools/NoahmpMacro.py` at
   configure time and requires Python 3; its source glob covers `src/`,
   `utility/`, and **`drivers/erf/` only** — `drivers/wrf/` is never compiled.
   New Fortran must land in `drivers/erf/` (or `src/`) to be built at all.

## 8. Verification

1. `make codegen-check` clean; `NoahmpIO_AssertAbi()` does not abort.
2. Build both ways — `cmake -DERF_ENABLE_NOAHMP=ON -DERF_ENABLE_NETCDF=ON` and
   `make USE_NOAHMP=TRUE USE_NETCDF=TRUE` — since both build systems list Noah-MP
   sources explicitly and must be updated together.
3. **Null test (the important one):** with `SF_URBAN_PHYSICS = 0` the existing
   `WPS_Test` case must reproduce its current trajectory **bitwise**. The urban
   code is entirely behind that guard.
4. **Degenerate-tile test:** with `SF_URBAN_PHYSICS = 1` and `FRC_URB2D ≡ 0`,
   results must match the `NATURAL_TABLE` rural run to round-off — this isolates
   the blend from the physics.
5. **Physics test:** urban case vs the Phase-0 WRF `sf_urban_physics = 1`
   reference. Compare the diurnal cycle of `TSK`, `HFX`, `LH`, and 2-m
   temperature over urban cells; expect the canonical urban signature — reduced
   latent flux, elevated storage, a nocturnal heat-island offset over the rural
   tile. Exact agreement with WRF is not expected (different dycore, different
   surface layer); a matching *sign and magnitude* of the urban–rural contrast is.
6. Unit tests from Phase 5 pass; the new regression case is added to
   `Tests/ERFGoldFiles`.
7. Restart round-trip: the urban prognostic state (`TR/TB/TG/TC`, `TRL/TBL/TGL`,
   `XXX*`) must be carried through checkpoint/restart via
   `NoahmpWriteRestartMod.F90` / `NoahmpReadRestartMod.F90`, or a restarted run
   will cold-start the canopy and diverge. Verify bitwise restart as
   `spec-noahmp-io.md` requires.
