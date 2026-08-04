# Plan: couple the Single-Layer Urban Canopy Model (SLUCM) into ERF

> Status: proposal, not yet implemented · Owns: the design and phasing for adding
> a WRF-style single-layer urban canopy model as an **urban tile** alongside the
> Noah-MP rural tile · Companions: [`spec-noahmp-api.md`](spec-noahmp-api.md)
> (the `NOAHMP` class and coupling fields), [`spec-noahmp-gpu.md`](spec-noahmp-gpu.md)
> (per-step staging), `Submodules/Noah-MP/drivers/erf/dev/spec-add-coupled-variable.md`
> (how a variable crosses the C++ ↔ Fortran boundary).
> Public docs: <https://erf.readthedocs.io/en/latest/CouplingToNoahMP.html>
>
> **Revision note.** This document was reviewed by three independent passes
> (fact-check, urban-physics, implementation-feasibility) and substantially
> rewritten. The single largest correction: the first draft asserted that SLUCM
> "blends in place into fields ERF already couples, so nothing new is required on
> the output side." **That is false in four ways** — see §1 and §4. Line and
> symbol references below were verified against ERF at this commit and the
> Noah-MP submodule at its pinned commit `e0aed20`.

## 1. Context — what we have, and what the tile blend actually is

ERF's real-data land surface is Noah-MP (`erf.land_surface_model = "NOAHMP"`,
selected at `Source/ERF.cpp:2469`, parsed at `Source/DataStructs/ERF_DataStruct.H:263`).
Over an urban cell it runs **bulk-urban**: `ConfigVarInTransferMod.F90:155-165`
sets `VegType = ISURBAN_TABLE` and `FlagUrban = .true.`, which swaps in an urban
roughness *and displacement height* (`src/GroundRoughnessPropertyMod.F90:70-76`)
and urban soil/thermal parameters. There is no canyon geometry, no separate
roof/wall/road energy balance, no urban thermal storage, no anthropogenic heat.

We want SLUCM (Kusaka et al. 2001; Kusaka and Kimura 2004; Chen et al. 2011 —
WRF `sf_urban_physics = 1`) run as a **tile** over the urban fraction, blended
with the Noah-MP rural tile by `FRC_URB2D`.

### What the SLUCM branch actually blends

Verbatim from the `IF (SF_URBAN_PHYSICS == 1)` branch of
`drivers/wrf/NoahmpUrbanDriverMainMod.F90` (lines 411-653) — **exactly nine
assignments**, at lines 565-581:

| Target | Form |
|---|---|
| `TS_URB2D` | `= TS_URB` (stored, not blended) |
| `ALBEDO` | `frc·ALB_URB + (1−frc)·ALBEDO` |
| `HFX` | `frc·SH_URB + (1−frc)·HFX` |
| `QFX` | `frc·LH_KINEMATIC_URB + (1−frc)·QFX` |
| `LH` | `frc·LH_URB + (1−frc)·LH` |
| `GRDFLX` | `frc·(G_URB·(−1.0)) + (1−frc)·GRDFLX` — **note the sign flip** |
| `TSK` | `frc·TS_URB + (1−frc)·TSK` |
| `QSFC` | `frc·QS_URB + (1−frc)·QSFC` |
| `UST` | `frc·UST_URB + (1−frc)·UST` |

**What is *not* blended, contrary to the first draft:**

- **`EMISS`** — never assigned in the SLUCM branch. The only `emiss(i,j) =` in the
  file is at line 796, inside the `SF_URBAN_PHYSICS == 2/3` block.
- **`TAU_EW` / `TAU_NS`** — the strings do not appear anywhere in the file (0 hits).
- **`ZNT`** — read in at line 447 (`ZNT_URB = ZNT(I,J)`) and passed as an input.
  `urban()` does set `ZNT = Z0` internally (`module_sf_urban.F:2400`), but the
  driver **never writes it back**.

The `*_RURAL` snapshot at lines 400-409 and the combined-radiation epilogue at
lines 744-831 (which recomputes `TSK` from `RL_UP_TOT` and `ALBEDO` from
`RS_ABS_TOT`) belong **exclusively to the BEP/BEM branch**. They are dead code
for SLUCM, which blends in place against the live `HFX(I,J)`. Do not copy them
as-is — but see §4, because ERF needs that epilogue's *physics* for a different
reason than WRF does.

## 2. Upstream survey — are we first?

Yes, for a mesoscale urban **parameterization**. From `erf-model/ERF` issues and
PRs (searched `urban`, `SLUCM`, `UCM`, `LCZ`, `canopy`, `building`, `land surface`):

- **No urban canopy parameterization exists and none is proposed.** No issue or
  PR mentions SLUCM, BEP, BEP+BEM, TEB, LCZ, or `FRC_URB2D`.
- **ERF's urban story upstream is building-resolving.** Immersed forcing plus
  embedded boundaries is actively developed: #2763, #2953, #3464/#3468, #3478,
  #3482 (open). `Docs/sphinx_doc/ERFvsWRF.rst:69` is the only `.rst` mention of
  "urban": *"Terrain and urban geometries may be simulated with immersed forcing
  or embedded (immersed) boundary techniques."*
- **Nearest precedent for a canopy parameterization**: the forest drag model
  (#1980, from AMR-Wind), `Source/SourceTerms/ERF_ForestDrag.{H,cpp}` applied at
  `ERF_MakeMomSources.cpp:747-783` — momentum sink only, no energy budget.
- **All 2026 LSM effort upstream is Noah-MP hardening**: #3414/#3412 (full
  coupling), #3255/#3317 (restart), #3286 (sea ice), #3334 (coarse→fine flux),
  #3428 (stress units), #3465 (2-D diagnostics).

Consequences: no branch to coordinate with, no reviewer with context — and ERF
already owns the *microscale* urban problem. SLUCM is the mesoscale complement,
and §7.1 must forbid running both over the same cells.

On the WRF side the tile driver now lives inside `NCAR/noahmp` as
`drivers/wrf/NoahmpUrbanDriverMainMod.F90`; only `phys/module_sf_urban.F` and
`run/URBPARM.TBL` remain in WRF proper.

## 3. What already exists in our tree

`Submodules/Noah-MP` pins `erf-model/noahmp` @ `e0aed20`. Much of the Fortran
urban plumbing is inherited from NCAR and never exercised (`SF_URBAN_PHYSICS = 0`):

| Present in `drivers/erf/` | What it gives us |
|---|---|
| `NoahmpIOVarType.F90-mc` (133 urban refs) | Most of the SLUCM array set — `tr/tb/tg/tc/qc/uc_urb2d`, `xxx{r,b,g,c}_urb2d`, `trl/tbl/tgl_urb3d`, `sh/lh/g/rn/ts_urb2d`, `psim/psih/u10/v10/GZ1OZ0/AKMS/th2/q2/ust_urb2d`, `cmcr/tgr/tgrl/smr/drel*/flxhum*`, `frc_urb2d`, `utype_urb2d`, `lp/lb/hgt/mh/stdh/lf_urb2d` — plus BEP/BEM and green-roof/PV sets, `sf_urban_physics`, `num_urban_*`, `urban_map_*`, `ISURBAN`, `URBTYPE_beg`, `LCZ_1..11_TABLE`, `IRI_URBAN`, `GMT`, `JULDAY`, `HRANG`, `DECLIN` |
| `NoahmpIOVarInitMod.F90-mc` (187 urban refs) | Allocation, gated on `SF_URBAN_PHYSICS == 1` vs `2/3`, plus a **second** guarded block holding the `undefined_real_neg` sentinel initialisation |
| `NoahmpReadNamelistMod.F90` | Reads `sf_urban_physics`, `use_wudapt_lcz`, `num_urban_hi`, `urban_atmosphere_thickness`, and validates the BEP/BEM level count |
| `NoahmpReadTableMod.F90` + `NoahmpTable.TBL` | `ISURBAN_TABLE`, `NATURAL_TABLE`, `LCZ_1..11_TABLE` (`URBTYPE_beg = 50`, LCZ 51-61) |
| **`ConfigVarInTransferMod.F90:155-165`** | **The rural-tile contract, already implemented.** `SF_URBAN_PHYSICS > 0` ⇒ `VegType = NATURAL_TABLE`, `GVFMAX = 96 %`, with the comment "urban is handled by explicit urban scheme outside Noah-MP" |
| `NoahmpInitMainMod.F90:175-207`, `NoahmpDriverMainMod.F90:178-194` | Urban-point LAI/biomass branching; urban irrigation |
| `NoahmpReadLandMod.F90` | Reads wrfinput directly (`XLAT`, `XLONG`, `IVGTYP`, `ISLTYP`, `VEGFRA`, `SHDMIN/MAX`) and the `ISURBAN` global attribute |
| `NoahmpIO.H-mc:169-170` | `Z0` and `ZNT` are **already coupled members**; Noah-MP fills `ZNT` (`EnergyVarOutTransferMod.F90:123`) |

**Corrections to the first draft's inventory:**

- The array set is **not complete**. Nine arrays that the WRF driver passes into
  `urban()` are absent from `NoahmpIOVarType.F90-mc`: `lf_urb2d_s` and the
  in-canyon vegetation set `tvg_urb2d`, `tt_urb2d`, `xxxvg_urb2d`, `tvgl_urb3d`,
  `smg_urb3d`, `cmcg_urb2d`, `flxhumvg_urb2d`, `flxhumt_urb2d`. In WRF these come
  from the registry-generated grid struct.
- `urban_map_*` is **declared and consumed as allocation extents but never
  assigned** at this commit (`grep "urban_map_[a-z]* *="` → zero hits). The
  namelist reader does *not* derive them.

### The four Fortran pieces that are missing

1. **`module_sf_urban.F90`** — the SLUCM physics, from WRF `phys/` (5620 lines).
   Must be renamed to `.F90`: both build systems glob `*.F90` only
   (`drivers/erf/CMakeLists.txt:39-44`, `Makefile:48`), and gfortran treats `.F`
   as **fixed-form** while the file is free-form. Its `use` list is
   `module_wrf_error` (line 7) and `module_model_constants, only: piconst`
   (line 11) — the file's `FATAL_ERROR`/`WRITE_MESSAGE` macros do *not* remove
   these. Brings `urban_param_init` (reads `URBPARM.TBL`) and `urban_var_init`.
2. **`NoahmpUrbanDriverMainMod.F90`** — from `drivers/wrf/`, reduced to the SLUCM
   branch. Its own `use` list needs `KARMAN/CP/XLV` rehoused and `cal_mon_day`
   resolved (§7.3). ~180 explicit dummy args — this is real plumbing, not a copy.
3. **`NoahmpUrbanShimMod.F90`** — supplying `wrf_message`, `wrf_error_fatal`,
   `piconst`, `KARMAN`, `CP`, `XLV`.
4. **`URBPARM.TBL`** (and `URBPARM_LCZ.TBL` if LCZ is enabled — a *separate*
   11-category file opened when `use_wudapt_lcz = 1`, `module_sf_urban.F:3018-3021`).

## 4. Gaps on the ERF (C++) side

The first draft claimed "nothing new is strictly required" on the output side.
**Wrong — but the gap is narrower and differently placed than the second draft
said.** `tau_ew`, `tau_ns`, `emiss` and all four banded albedos are **already
coupled end to end** (`ERF_NOAHMP_Fields.H:49-53,77-81`; `NoahmpIO.H-mc:41-47`;
written by `EnergyVarOutTransferMod.F90:107-110`). ERF receives them from
Noah-MP today. What is missing is that **SLUCM never contributes to them** — so
every fix below is *Fortran-side blend work in Phase 2*, not ERF-side coupling
work, and none of it touches the ABI or `m_lsm_data_size`.

The single genuine coupling gap is **`QFX`**: 0 hits in `NoahmpIO.H-mc`,
Fortran-only at `NoahmpIOVarType.F90:315`. It must be promoted in Phase 3
(appending to `NOAHMP_OUTPUT_2D_FIELDS_TAIL` does **not** change
`m_lsm_data_size`, so no checkpoint break).

| Gap | Detail |
|---|---|
| **Momentum never reaches ERF** | WRF conveys urban drag through `UST` because its surface-layer scheme reconstructs stress from `u*`. **ERF does not**: `SurfaceLayer` reads `tau13`/`tau23` from the LSM directly (`ERF_SurfaceLayer.cpp:658-735`) and back-derives `u* = sqrt(|τ|)` (`:1063`); `UST` is not a coupled field. Ported literally, the atmosphere sees **rural stress at full weight** with an **urban-weighted heat flux**. Over a `frc = 0.9` city (`Z0C ≈ 0.5 m` vs rural `z0m ≈ 0.06-0.1 m`) τ is under-predicted ~1.8×, and `t* = −t_flux/u*` then pairs an urban `t_flux` with a rural `u*` — over-estimating `t*` ~2-3× and handing the PBL a spurious super-unstable surface layer over every city at midday. **Must add a `TAU_EW/NS` blend** from `UST_URB² = (R·CDR + RW·CDC)·UA²` (kinematic — multiply by `RHOO`), projected on `U1_URB/UA_URB`, `V1_URB/UA_URB`. |
| **Urban albedo never reaches radiation** | RRTMGP consumes four *banded* albedos (`sfc_alb_dir_vis/nir`, `sfc_alb_dif_vis/nir`, `ERF_Radiation.cpp:620-622`), filled from `ALBSFCDIRXY`/`ALBSFCDIFXY`. SLUCM's `ALB_URB` is blended only into the scalar `ALBEDO`, which in ERF is the **diagnostic** `o_albedo`. As written, an urban cell reflects shortwave with the *rural* albedo — the primary UHI forcing is silently dropped. **Must blend all four bands** (SLUCM is broadband; apply `ALB_URB` to all four and document), and hold `alb_rural` when `swdown ≤ 0` since `urban()` zeroes `ALB` at night (`module_sf_urban.F:2384`). |
| **Emissivity and the meaning of `TSK`** | `EMISS` is not blended at all, and `TS_URB` is **not a radiative temperature** — it is `TA + FLXTH/CHS` (`module_sf_urban.F:2446`), a bulk-flux inversion; the radiative form is commented out one line later. ERF's rural `TSK` *is* radiative (`TemperatureRadSfc`). So `TSK_blend` mixes two different quantities and is handed to RRTMGP with a rural `sfc_emis`, while SLUCM's true upward longwave `LW_URB` (including canyon trapping) is discarded. The T⁴ linearisation error is small (~0.2-0.4 K); the emissivity/flux-temperature mismatch is tens of W m⁻². **Must port the BEP epilogue's radiation algebra to the SLUCM branch** — `rl_up_urb = −LW_URB`, blend `rl_up` and `emiss`, invert for a genuine radiative `TSK`; keep `TS_URB2D` as the WRF-comparable diagnostic. |
| **Latent-heat constant mismatch** | SLUCM uses `ELL = 2.442e6`; ERF recovers `q_flux = LH/(ρ·L_v)` with `L_v = 2.5e6` (`ERF_Constants.H:59`) — the urban moisture flux arrives **2.3 % low** (Noah-MP's rural tile is 0.4 % high in the other direction). **Couple `QFX` directly** and take `q_flux = QFX/ρ`; the SLUCM branch already blends `QFX`. Fixes both tiles. |
| **Roughness** | `ZNT` already crosses the ABI and Noah-MP fills it — so wiring `ZNT` → a new `LsmData` row → `RoughCalcType::LSM` is genuinely ERF-only and unblocked. **But that carries the *rural/bulk* roughness only**, because the SLUCM branch never writes `ZNT` back. Adding the urban blend requires a submodule change, and a *linear* blend of `z0` is wrong — the flux-consistent form averages `1/ln²(z/z0)` at blending height. Meanwhile `rough_type_land` hard-`Abort`s unless `"constant"` (`ERF_SurfaceLayer.H:211-215`), and ERF has **no displacement height** (`ERF_MOSTStress.H:174` is `log(zref/z0)`) while `ZDC ≈ 0.76·ZR` is 4-8 m for the default table. |
| **Static urban input** | `FRC_URB2D`, `URB_PARAM`, `LU_INDEX`, LCZ read nowhere in ERF. **But this is optional**: `urban_var_init` derives `UTYPE_URB2D` from `IVGTYP` and falls back to `FRC_URB_TBL(UTYPE)` when `FRC_URB2D` is absent, and to table geometry when `HGT_URB2D ≤ 0`. Phase 1 can run without new static data. |
| **`urb_frac_lev`** | Declared (`ERF.H:1033`), allocated `setVal(one)` (`ERF_MakeNewArrays.cpp:501-504`) — **never filled, never read** (repo-wide grep returns only those 6 lines). Note the default is **1.0 = fully urban**, the wrong default the moment anything reads it. |
| **Solar geometry** | `orbital_decl` supplies declination (`ERF_OrbCosZenith.H:20`). **There is no hour-angle output** — `h` in `orbital_avg_cos_zenith` is the half-day length. `OMG_URB2D` is new code: `OMG = 2π·frac(jday) + lon_rad − π`, verified against `tloc` giving 12 at local solar noon. `XLAT` must reach `urban()` in **degrees**. |
| **Advance dispatch** | `ERF_AdvanceLSM.cpp:13` hard-codes `LandSurfaceType::NOAHMP`. **Non-issue under Design A** — noted only to stop someone adding a `LandSurfaceType::SLUCM`. |

## 5. Architecture decision

**Recommendation: SLUCM as an urban tile inside the Noah-MP Fortran driver**
(design A), mirroring `drivers/wrf/NoahmpUrbanDriverMainMod.F90`. Reject an
ERF-native C++ scheme (design B) because:

- The rural-tile contract already exists in `ConfigVarInTransferMod.F90:155-165`.
- `LandSurface` holds **one model per level** and `SetModel<T>()`
  (`ERF_LandSurface.H:29-36`) assigns the same type to every level, so
  "Noah-MP for the rural tile, SLUCM for the urban tile" is unreachable by
  construction. *(The first draft cited the `typeid` guard at
  `ERF_LandSurface.H:153-165`. That guard compares two
  `std::unique_ptr<NullSurf>` — non-polymorphic, resolved statically, always
  true. The conclusion holds; the mechanism cited did not.)*
- `TSK`/`ALBEDO`/`EMISS` feed RRTMGP; blending in Fortran keeps one source of truth.
- Bit-level traceability to WRF is the validation strategy (§8.4).

Costs accepted: vendoring WRF Fortran into the fork, and SLUCM staying host-side
(matching Noah-MP, which already round-trips device→host each step).

**Where design B would be right**: a multi-layer scheme (BEP/BEP+BEM) injecting
drag and TKE at several levels is not an LSM tile — it belongs in
`ERF_MakeMomSources.cpp` like `ForestDrag`. Future work; the `a_u_bep`,
`a_t_bep`, `sf_bep`, `vl_bep` arrays are already declared.

## 6. Phased plan

Reordered from the first draft: unblocked, independently-useful ERF-side work and
the missing test infrastructure now come **first**, because Phases with "Fortran"
in them are all inside `erf-model/noahmp`.

### Phase 0 — decisions and prerequisites (no SLUCM code)

Four decisions, all of which change later code, plus two prerequisite bug fixes:

1. **`CHS` is settled — do not defer.** `CHS` is `INTENT(INOUT)` but never
   modified inside `urban()`; it is used only in `TS = TA + FLXTH/CHS` and
   `QS = QA + FLXHUM/CHS`. `SH`, `LH`, `G`, `RN`, `UST` and all prognostic state
   are independent of it. Pass `NoahmpIO%CHXY` (Noah-MP's `ExchCoeffShSfc` —
   same quantity, same units `[m/s]`, same column, same reference height),
   floored at `1.0e-2` as WRF does. `CHS2 = FVEG·CHV2XY + (1−FVEG)·CHB2XY`.
   Do **not** enable the commented-out internal `CHS` lines — that breaks
   bit-comparability with WRF. Combined with the radiative-`TSK` fix (§4), `CHS`
   stops affecting anything ERF consumes for radiation.
   *(Note: `SFCDIF` is not a `URBPARM.TBL` key. The real knobs are `CH_SCHEME`
   (1 = M-O via `mos`; 2 = Narita, the default) and `AKANDA_URBAN`. `SFCDIF_URB`
   is always used above canopy; `louis79`/`louis82` are dead code.)*
2. **Reference height.** Decide that ERF stages `ZLVL` from the true `z(klo)`
   (§7.2) rather than the namelist constant, and fix the resulting
   `RefHeightAboveSfc` inconsistency in Noah-MP.
3. **Precision strategy** (§7.5) — per-source promotion flags, not a 5620-line
   kind rewrite.
4. **Supported envelope** (§7.2) — `dx ≳ 1 km` and `z₁ ≥ 2·max(ZR)`, aborting outside.

Prerequisite fixes, each its own PR:

- **Wire `JULIAN`/`YR`/`GMT`.** `NoahmpIOVarInitMod.F90-mc:882-883` hard-sets
  `YR = 2000`, `JULIAN = 1.0` ("wire via ABI for correctness") and nothing ever
  updates them. This already breaks Noah-MP phenology under
  `DYNAMIC_VEG_OPTION = 4` (what `WPS_Test` uses) and silently disables the
  existing urban irrigation. SLUCM's entire shadow/solar model keys off it.
- **Hoist `CAL_MON_DAY`** from `NoahmpDriverMainMod.F90:224` into `utility/`
  (§7.3).

### Phase 1 — ERF-side work that needs no submodule change

Deliverable in `sunt05/ERF` alone, and useful with or without SLUCM:

1. **CI that builds Noah-MP at all.** `grep -i noahmp .github/workflows/` returns
   nothing across all 18 workflows, so `Tests/Unit/LandSurfaceModel/Noah-MP/ERF_GTestNoahMPResultPolicy.cpp`
   — behind `if(ERF_ENABLE_NOAHMP)` at `Tests/Unit/CMakeLists.txt:149` — **has
   never compiled in CI**. Clone `gcc-rrtmgp.yml`, add `gfortran`,
   `libnetcdff-dev`, `python3`, `submodules: recursive`, `-DERF_ENABLE_NOAHMP=ON`.
   Everything else depends on this existing.
2. **Regression baseline — blocked on a data decision, not on code.**
   `Tests/CTestList.cmake` has no Noah-MP or WPS entry. (`WPS_Test` needs **no**
   `GNUmakefile`: ERF builds a single shared `erf_exec`, and every `add_test_r`
   entry passes `TEST_DIR=""`.) The real blocker: `wrfinput_chisholmview_d01`,
   `wrfbdy_...` and the four RRTMGP coefficient files are **not in the repo**;
   gold files are 334 MB committed in-tree with **no download mechanism in any
   GitHub workflow** (the gold-file repo is LLNL-GitLab-only). The nearest
   precedent, `add_test_r(Radiation)`, is itself dead behind an
   `ERF_ENABLE_RRGMTP` typo at `CTestList.cmake:523`. **Decide where test data
   lives, and size a runner-appropriate deck, before writing this.**
3. **`ZNT` → MOST roughness.** Append `X(znt)` / `X(o_znt, ZNT)` /
   `X(ZNT_o, znt, o_znt)` to the three registries in `ERF_NOAHMP_Fields.H` —
   `read_results` is fully table-driven, so **`ERF_NOAHMP_Advance.cpp` needs no
   edit**. Add a positivity rule for `o_znt` in `ERF_NOAHMP_ResultPolicy.H`
   (`result_is_valid` accepts `0.0`, and `ERF_MOSTStress.H:174` computes
   `log(zref/z0)` → `inf`). Add `RoughCalcType::LSM`, accept
   `erf.most.roughness_type_land = "lsm"` (`ERF_SurfaceLayer.H:211-215`), and
   **widen the six `Abort("Unknown value for rough_type_land")` guards** rather
   than writing new functors. `z_0[lev]` has no per-step land update path today
   (filled once at init, updated only for `!is_land`) — add `get_lsm_z0()`
   mirroring `get_lsm_tsurf` (`ERF_SurfaceLayer.cpp:1214`) plus an
   `m_lsm_z0_indx` resolved by name. **Extract the per-cell selection into a
   header-only free function** (as `ERF_SurfaceLayerStress.H` does) — otherwise
   it is not unit-testable.
   Two caveats: this **breaks existing checkpoints in Phase 1**, because
   `NOAHMP_LSMDATA_FIELDS` *is* `m_lsm_data_size` and LSM data is checkpointed
   positionally (`ERF_Checkpoint.cpp:238-244`, `:818-825`) with the soil profile
   indexed off `NumVars`. And it is **largely inert where SLUCM runs**: over land
   with valid LSM fluxes on both neighbouring cells, ERF discards MOST entirely
   (`ERF_SurfaceLayerStress.H:62-64`) and re-derives `u* = sqrt(|τ|)`. Its value
   is at LSM-invalid and land/sea-boundary faces, and for PBL `get_z0` consumers.
4. **`urb_frac_lev` becomes live.** Add `FRC_URB2D` to
   `ERF_ReadFromWRFInput.cpp:125-129` and to `NC_names`. **Critically, also add
   it to the `has_fallback_behavior` disjunction at
   `ERF_InitFromWRFInput.cpp:352-355`** — that list is nine dycore variables, and
   anything else failing to read hits `amrex::Abort`, so without this **every
   existing wrfinput (none of which carry `FRC_URB2D`) dies at init.** Fill
   `urb_frac_lev` with a clamped copy modelled on the `TSK` block, change the
   default from `setVal(one)` to zero (`ERF_MakeNewArrays.cpp:503`), and register
   it in `ERF_Plotfile2DCatalog.cpp` (Geometry category — needs an explicit fill
   block in `ERF_Plotfile2D.cpp` in **catalog order**, immediately after
   `landmask`). **Drop `LU_INDEX`** — redundant with `IVGTYP`, no ERF container,
   pure added abort risk.
5. **Parameter-check aborts** for SLUCM with `ImmersedForcing` (both
   `BuildingsType` and `TerrainType`) or `ForestDrag`. Home is
   `ERF::ParameterSanityChecks()` (`ERF.cpp:2492`), not `SolverChoice::init_params`,
   because `do_forest_drag` is only set later at `ERF.cpp:2442`; the existing
   idiom is the `cf_width` check at `:2551`. Precedent: ERF already errors on
   ForestDrag + ImmersedForcing at `ERF_MakeMomSources.cpp:155-156`.
   **Blocked in Phase 1**: `sf_urban_physics` lives only in the Fortran namelist
   and is not an ABI member, so ERF cannot know SLUCM is on. Either add an
   `erf.slucm` mirror flag here, or defer this item to Phase 3 and check after
   `NOAHMP::Init`.
6. **Delete dead build-graph forcing.** *(The second draft had this backwards.)*
   `Exec/GNUmakefile:48-52` includes `Make.ERF.general` **only when
   `USE_RRTMGP=TRUE`**, so its `USE_NOAHMP ⇒ USE_RRTMGP` override at
   `Make.ERF.general:36-38` is unreachable for the case it claims to guard —
   `make USE_NOAHMP=TRUE` goes through `Make.ERF`, which requires only NetCDF,
   exactly matching CMake. CMake and `Make.ERF` already agree; delete the dead
   lines rather than mirroring them.
7. **Unit tests** for 3 and 4, registered in the now-live `if(ERF_ENABLE_NOAHMP)`
   block. Testable: the `ResultPolicy` `o_znt` rule (the existing X-macro loop in
   `ERF_GTestNoahMPResultPolicy.cpp` covers the new row for free), the extracted
   roughness-selection free function, and the plotfile catalog entry plus a
   catalog↔fill-order regression. **Not** testable in isolation: `SurfaceLayer`
   itself (no test in the tree constructs one) and the NetCDF read path.

### Phase 2 — Fortran: vendor SLUCM into the submodule fork

Branch off `e0aed20`. **Both build systems glob**, so no build-file edits are
needed for new `*.F90` in `drivers/erf/` — but two naming constraints bind:
the extension must be `.F90`, and **the basename must equal the module name**,
because the Makefile's dependency generator emits `obj/<module-in-use-stmt>.o`
for every `use` (`Makefile:83-84`).

1. `drivers/erf/module_sf_urban.F90` — renamed, `use` lines redirected to the
   shim, with a small **documented patch list** (§7.6), each marked `! ERF: <reason>`.
2. `drivers/erf/NoahmpUrbanShimMod.F90` — `wrf_message`, `wrf_error_fatal`,
   `piconst`, `KARMAN`, `CP`, `XLV`. **Do not** name it `module_model_constants`:
   ERF already ships `Source/Microphysics/Morrison/ERF_module_model_constants.F90`
   declaring that module with its own `piconst`, and the GNU make path puts the
   submodule's `include/` (which receives every `.mod`) on ERF's
   `INCLUDE_LOCATIONS` — so `USE_MORR_FORT=TRUE` + `USE_NOAHMP=TRUE` can resolve
   the wrong `.mod`.
3. `drivers/erf/NoahmpUrbanDriverMainMod.F90` — SLUCM branch only, BEP/BEM
   branches and their `use` lines dropped. **Add** the outputs from §4:
   `TAU_EW/NS` blend, four banded albedos, `EMISS` + an urban effective
   emissivity, radiative `TSK` via the ported radiation algebra. Snapshot
   `TAU_EW/NS` and `EMISS` into the `*_RURAL` block. Delete the hard-coded
   `IF (I.EQ.73.AND.J.EQ.125)` debug point at lines 530-532.
4. Declare the nine missing arrays (§3) if `TREEOPTION`/distributed aerodynamics
   are ever enabled; otherwise scope them out and assert the options are off.
5. Call the driver from `NoahmpDriverMainMod.F90` under
   `if (SF_URBAN_PHYSICS == 1)`. Call both initialisers from
   `NoahmpInitMainMod.F90`, but **latch only `urban_param_init`** — it `OPEN`s
   `URBPARM.TBL` and fills module-level `SAVE`d `*_TBL` arrays, so once per
   process is right, and `InitMain()` runs per box (`ERF_NOAHMP_Init.cpp:195`).
   `urban_var_init` is a **per-tile** initialiser (`DO I=ims,ime / DO J=jms,jme`
   writing `UTYPE_URB2D`, the `FRC_URB2D` fallback, `TR/TB/TG/TC`,
   `TRL/TBL/TGL_URB3D`, `module_sf_urban.F:3711-3900`) — latching it would leave
   every box after the first with uninitialised urban state.
9. **Urban restart I/O.** `NoahmpWriteRestartMod.F90` / `NoahmpReadRestartMod.F90`
   are hand-enumerated varid lists with **zero** urban entries. Adding the §8.10
   state set is a named Phase-2 deliverable, and it is not small.
6. Guard `SF_URBAN_PHYSICS == 1 .and. NSOIL /= 4` — `urban_param_init` forces
   `num_roof/wall/road_layers = num_soil_layers` and `urban_var_init` hard-codes
   `TRL_URB3D(I,1..4,J)`.
7. Set `VEGFRA = 0.96*100` alongside `GVFMAX` in `ConfigVarInTransferMod.F90`
   (§7.4), or abort for `OptDynamicVeg ∈ {1,6,7}`.
8. Ship `URBPARM.TBL`; ship `URBPARM_LCZ.TBL` or gate `USE_WUDAPT_LCZ = 0`.

### Phase 3 — interface (coupled variables)

Per `spec-add-coupled-variable.md`, with four corrections the first draft missed:

- **Always append** to the `@NoahmpMacro:Source` block. Order *is* ABI order, the
  mirror is a flat pointer array, and `NoahmpIO_AssertAbi()` checks **only
  element precision** — there is no member-count check. With a half-stale build
  (the GNU path links a prebuilt `libnoahmp.a` and `cp -u`s the header), an
  appended member is merely garbage while a **middle insertion silently corrupts
  every later slot**. Add `NoahmpIO_NumMembers_fi()` to the ABI assert.
- **`utype_urb2d` cannot cross today.** `tools/NoahmpMacro.py` matches only
  `NoahmpArray[23]D<noahmp_real>` and `KIND_TRAITS["array"]` hardcodes
  `real(kind=c_kind_noahmp)` — there is no integer-array kind (which is why
  `IVGTYP`/`ISLTYP` are hand-written). Either extend the generator with an
  `iarray` kind (one regex + one traits row; `NoahmpArray.H` is already
  `template <typename T>`) or keep it Fortran-side.
- **Promotion turns a guarded allocate into an unconditional one.** The urban
  arrays are allocated inside `if (SF_URBAN_PHYSICS > 0)` and their sentinel
  initialisation sits in a **second** guarded block. Promote and the generated
  allocate becomes unconditional while the sentinel init stays guarded — leaving
  them **allocated but uninitialised** in every non-urban run, so C++ reading
  `FRC_URB2D` gets garbage rather than a detectable sentinel. Move the sentinel
  assignments to the unconditional block, or do not promote.
- **Scalars have no allocate to delete.** `JULIAN`/`GMT`/`DECLIN`/`YR` become
  C++-owned pointers on promotion, so every existing Fortran write becomes a
  write through that pointer and must be sequenced against `ScalarInitDefault()`.

**Promote `QFX`** (§4): append to the `@NoahmpMacro:Source` block and add
`X(o_qfx, QFX)` to `NOAHMP_OUTPUT_2D_FIELDS_TAIL`. This is the only genuinely
missing coupled output; it does not change `m_lsm_data_size`, so no checkpoint
break. If Phase 1 item 5 was deferred, promote `sf_urban_physics` as an int
scalar here too, and move the cross-option abort to after `NOAHMP::Init`.

Settle the static-field design **once**: static urban fields stay Fortran-side
via `ReadLandMain`; only `FRC_URB2D` crosses, once, out, at init. `XLONG` need
not cross — `NoahmpReadLandMod.F90:187` already reads it. `NOAHMP_INPUT_2D_FIELDS`
is the wrong home regardless: it is the **per-step** staging table, so static
morphology would be re-copied every timestep.

### Phase 4 — ERF driver

- `ERF_NOAHMP_Fields.H`: only `QFX` needs a new row — `tau_ew`, `tau_ns`,
  `emiss` and the banded albedos are already registered (§4). **The
  checkpoint-invalidating change is `X(znt)` in Phase 1**, not here, because
  `NOAHMP_LSMDATA_FIELDS` is what sizes `m_lsm_data_size`; `QFX` goes only into
  the output table. Add a stored `m_lsm_data_size` to the checkpoint header with
  a mismatch abort, or document the break.
- `ERF_NOAHMP_Advance.cpp`: stage `DECLIN`, `OMG`, `JULIAN`, `YR`, `GMT`, and the
  true `ZLVL`; take `q_flux = QFX/ρ`.
- `ERF_NOAHMP_Init.cpp`: copy `FRC_URB2D` into `urb_frac_lev`; require a per-level
  land file when SLUCM is on (urban morphology must not be interpolated across
  levels by `interp_from_lev0`).

### Phase 5 — validation

See §8.

## 7. Risks and open questions

1. **Double counting.** Abort if SLUCM is active with `buildings_type = ImmersedForcing`
   or `ForestDrag` on the same level.
2. **Reference height — the real problem is inconsistency, not the guard.**
   `ZA_URB = 0.5·DZ8W` and `DZ8W = 2·ZLVL`, so `ZA_URB ≡ ZLVL` — the **namelist
   constant** (`WPS_Test` uses `ZLVL = 10.0`), not where ERF samples forcing
   (`k = klo`, cell centre **46.875 m** in that case). A third value,
   `erf.most.zref = 1.0`, sits in the inputs file. With the shipped `URBPARM.TBL`
   (`ZDC ≈ 0.762·ZR`), `urban()` hard-`FATAL_ERROR`s when `ZA ≤ ZDC + Z0C + 2`:

   | urban type | ZR | ZDC | Z0C | aborts if `ZA ≤` | canyon-wind fallback if `ZA ≤` |
   |---|---|---|---|---|---|
   | 1 low-density residential | 5.0 | 3.81 | 0.172 | 5.98 m | 7.0 m |
   | 2 high-density residential | 7.5 | 5.72 | 0.333 | 8.05 m | 9.5 m |
   | 3 commercial | 10.0 | 7.62 | 0.531 | 10.16 m | 12.0 m |

   **At the configuration Phase 1 ships, none of these fire.** With
   `use_wudapt_lcz = 0`, `urban_var_init` assigns `UTYPE_URB2D = 2` to *every*
   `IVGTYP == ISURBAN` cell (`module_sf_urban.F:3745-3749`); types 1 and 3 are
   reachable only through LCZ vegtypes, which §4 defers. So `ZLVL = 10.0` clears
   type 2's 8.05 m abort *and* its 9.5 m fallback. **The abort and the fallback
   go live only with LCZ/NUDAPT morphology or an `MH_URB` override.** Justify
   staging `ZLVL` on reference-height consistency, not on a phantom crash.
   Fix: stage `ZLVL` from `z(klo)`; require `ZA ≥ max(ZDC + Z0C + 2, ZR + 2)`
   and warn unless `ZA ≥ 2·ZR` (MOST is valid only above the roughness sublayer,
   2-5 building heights). **Supported envelope: `dx ≳ 1 km` and `z₁ ≥ 2·max(ZR)`;
   abort outside it.** At grey-zone `dx` users refine `Δz` to O(10-20 m), putting
   `z₁` inside the canopy, and at `dx ≲ 500 m` the horizontal-homogeneity
   assumption behind `FRC_URB2D` fails anyway.
3. **Circular module dependency.** The urban driver needs `cal_mon_day`, which
   lives inside `NoahmpDriverMainMod` — the module that will *call* the urban
   driver. `use NoahmpDriverMainMod` there is a hard gfortran error and a literal
   cycle in the generated `.d` files. Hoist it to `utility/` first.
4. **Rural tile is under-specified.** `GVFMAX = 96 %` only works for
   `OptDynamicVeg ∈ {2,3,4,5,8,9}`; for `{1,6,7}` the rural tile inherits the
   urban cell's greenness (5-20 %) and the mosaic is meaningless. Also: LAI from
   the wrfinput time series is near zero over urban cells under
   `DYNAMIC_VEG_OPTION ∈ {7,8,9}`; `ISLTYP` is often a placeholder; **SLUCM has
   no snow at all** (no roof/road snow albedo or melt) so winter cases must not
   be used for Phase-1 validation; `RAIN_URB` treats snow as liquid; and
   `SFCRUNOFF`/`UDRUNOFF`/`SMSTOT` are whole-cell values that physically exist
   only over `(1−frc)`.
5. **Precision.** The submodule builds `DOUBLE_PREC`/`c_kind_noahmp`;
   `module_sf_urban` uses bare `REAL`. A 5620-line kind rewrite contradicts
   keeping the vendored diff readable. Use per-source promotion flags —
   gfortran `-fdefault-real-8 -fdefault-double-8` (**both**, or `DOUBLE PRECISION`
   becomes 16 bytes), Intel/NVHPC `-r8`, Cray `-s real64` — via
   `set_source_files_properties` and a file-specific make rule. Guard the
   single-precision build.
6. **Latent defects in the vendored revision.** This revision of
   `module_sf_urban.F` has uninitialised reads (`TGEP = TGE`, `SROOTP = SROOT`;
   `IRI_SCHEME = 1` reading `tloc` assigned only under `ahoption == 1`; `SW`
   reused as street width), and a `print*` inside the `TREEOPTION = 1` shortwave
   branch that emits one line per urban column per step. All are dormant at the
   shipped defaults but argue for a documented patch list rather than "unmodified".
7. **Default option set.** `AHOPTION`, `ALHOPTION`, `GROPTION`, `TREEOPTION`,
   `IRI_SCHEME` are **all 0** by default and `IMP_SCHEME = 1`. Phase 1 ships
   exactly those — the WRF-comparable configuration, and the only one whose code
   paths avoid the §7.6 defects. Expose `AHOPTION` (anthropogenic heat) as a
   documented knob since it needs only `AH_TBL`/`AHDIUPRF` plus a correct `OMG`;
   keep `TREEOPTION`/`GROPTION` out until upstream fixes the uninitialised state.
8. **Direct/diffuse partitioning is hard-coded.** `module_sf_urban.F:948-951`
   sets `SSGD = 0.75·SSG` unconditionally, discarding both the driver's argument
   and ERF's real banded direct/diffuse fluxes from RRTMGP. Under overcast skies
   SLUCM casts shadows that do not exist. Route ERF's real direct fraction in
   behind a flag defaulting to WRF behaviour.
9. **`urban_var_init` on restart clobbers state.** `QC_URB2D = 0.01` sits
   **outside** the `IF (.not.restart)` guard. Decide explicitly whether ERF calls
   it on restart; if so, move that line inside the guard.
10. **Thread safety.** `module_sf_urban`'s module-level `SAVE` state is not
    thread-safe. Safe today (`ERF_NOAHMP_Advance.cpp:387`'s `MFIter` loop has no
    OpenMP pragma, `USE_OMP = FALSE`), but it blocks future OpenMP tiling.
11. **Licensing.** WRF is public domain, so vendoring in is clean; but the
    combined fork stays under the UCAR Noah-MP license, whose no-resale clause is
    a redistribution constraint on BSD-3 ERF. Status quo (Noah-MP is already a
    submodule) — but get it confirmed rather than asserted.
12. **Vendoring drift.** Record the WRF source revision in the file header so a
    future re-sync is a readable diff.
13. **Submodule ownership.** Phases 2-3 are entirely inside `erf-model/noahmp`.
    Work on a personal fork (`.gitmodules` `url` is a one-line change and
    `shallow = false` already), then PR upstream and flip the URL back before merge.

## 8. Verification

`make codegen-check` is **not** a verification: the Makefile runs the generator
at parse time (`$(info … $(shell python3 …))`) and CMake at configure time, so
the check always compares against files it just rewrote, and the targets are
gitignored. Either drop it or make the Makefile skip the parse-time run when
`MAKECMDGOALS` is `codegen-check`.

1. **Build both paths** — CMake and GNU make — and confirm `NoahmpIO_AssertAbi()`
   does not abort.
2. **Null test.** With `SF_URBAN_PHYSICS = 0`, `WPS_Test` reproduces its
   trajectory **bitwise**. Requires the Phase-1 baseline to exist first.
3. **Per-column energy-balance closure** (cheapest, highest value). `urban()`
   defines `G = −FLXG·697.7·60` (`module_sf_urban.F:2386`) and
   `RN = (SNET+LNET)·697.7·60` (`:2387`), while the per-facet residual is
   `SR + RR − HR − ELER − G0R = 0`. The closing identity is therefore
   **`RN = SH + LH − G`**, i.e. assert
   `|RN_URB − SH_URB − LH_URB + G_URB| < 1e-6·max(1,|RN_URB|)`.
   *(The second draft wrote `− G_URB`, which is off by `2G` — O(100 W m⁻²)
   mid-afternoon — and would fail on a **correct** port.)*
   Gate the assert on `AHOPTION == 0 .and. ALHOPTION == 0`: anthropogenic heat
   enters `FLXTH`/`FLXHUM` but not `RN`, so with §7.7's advertised `AHOPTION`
   knob the identity becomes `RN = SH + LH − G − AH − ALH`.
4. **Grid-cell closure across the blend.** `LW_URB = LLG − LNET·697.7·60` is the
   **total** upward longwave (emission *plus* reflection), so the net is
   `GLW − LW_up`, not `ε·GLW − LW_up`. Assert
   `SWDOWN·(1−α_blend) + GLW − LW_up_blend = HFX + LH + GRDFLX`.
   **This fails against the first draft's design** — it is precisely what exposes
   the missing banded albedo and emissivity, and why §4's radiation fix is
   mandatory rather than optional.
5. **Offline single-column benchmark.** Drive the vendored `urban()` with the
   Grimmond et al. (2010, 2011) international urban energy-balance comparison
   forcing (Vancouver–Sunset, Melbourne–Preston) and compare `Q*`, `Q_H`, `Q_E`,
   `ΔQ_S` against the published SLUCM submission and the observations. This is
   the only test that separates "we broke the port" from "we broke the coupling".
6. **Single-column WRF↔ERF bit-comparison.** Identical forcing arrays, identical
   `URBPARM.TBL`, one column, 48 h; require `SH/LH/G/TS/TR/TB/TG` agreement to
   `< 1e-10` relative in double precision. Achievable because it is the same
   Fortran, and it is the real regression guard on the vendoring and the kind
   handling. Keep a 3-D WRF comparison as a sanity check, **not** as the
   acceptance gate — dycore differences swamp it.
7. **Momentum coupling** — nothing else in this list tests §4's first blocker.
   Assert `τ_blend = frc·ρ·UST_URB² + (1−frc)·τ_rural`, and that at `frc = 1`
   the `u*` ERF recovers via `sqrt(|τ|)` equals `UST_URB`.
8. **`frc` sweep.** `frc ∈ {0, 0.25, 0.5, 0.75, 1}` on a single column; assert
   every *blended* quantity is linear in `frc` and that `frc = 1` reproduces the
   pure-urban run. Note `u*` is **not** linear in `frc` by construction
   (`u* = sqrt(|τ|)`, `ERF_SurfaceLayer.cpp:1063`) — assert linearity of `τ`. *(The first draft's "set `FRC_URB2D ≡ 0`" test
   is not executable: `urban_var_init` detects `FRC_URB2D ≤ 0` on an urban
   `IVGTYP` and overwrites it with `FRC_URB_TBL(UTYPE)` = 0.5/0.9/0.95. To test
   the null path, set every urban `IVGTYP` to `NATURAL_TABLE` instead.)*
8. **Reference-height sensitivity sweep.** `ZA ∈ {6,10,15,20,30,47,60}` m per
   urban type, tabulating `ΔSH`, `ΔTS`, `Δu*`. Turns §7.2's envelope from a guess
   into a documented number.
9. **Spin-up.** `TRL/TBL/TGL` initialise from `TSK`/`TSLB` — a poor guess for
   concrete at `CAPR = 1e6 J m⁻³ K⁻¹`. Require day-N → day-N+1 drift in `TBL(4)`
   below 0.1 K and state a minimum spin-up (5-10 days) before any observational
   comparison, or §8.5's "nocturnal heat-island offset" measures initialisation.
10. **Restart, per variable.** Dump and diff **every** state variable across the
    restart boundary, not just the trajectory. The minimal bitwise set:
    `tr/tb/tg/tc/qc_urb2d`; `trl/tbl/tgl_urb3d(1:4)`; `xxx{r,b,g,c}_urb2d`;
    **`cmr/chr/cmc/chc_sfcdif`** (the `SFCDIF_URB` iteration state, relaxed with
    `WOLD = 0.15` — genuinely prognostic and missed by the first draft);
    `flxhum{r,b,g}_urb2d` and `drel{r,b,g}_urb2d` (their `*P = *` prologue is
    unconditional). Add `cmcr/tgr/tgrl/smr` if `GROPTION = 1`, and the tree set
    if `TREEOPTION = 1`. `uc_urb2d` is diagnostic (recomputed from `UA` before
    first use). None of these appear in `NoahmpWriteRestartMod.F90` today.
