/*
 * NOAHMP::Init: builds the surface (lsm) geometry and coupling MultiFabs, sizes
 * and initializes one NoahmpIO_type per box, and broadcasts the firing parameters.
 */

#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <ctime>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>

#include <ERF_NOAHMP.H>
#include <ERF_Constants.H>
#include <ERF_EpochTime.H>
#include <NoahmpFatal.H>

using namespace amrex;

#ifdef ERF_USE_NETCDF
// Defined in Source/Initialization/ERF_InitFromWRFInput.cpp. Declared here the
// same way ERF.cpp:130 and ERF_Tagging.cpp:8 declare it -- it has no header.
// Collective: reads on the IOProcessor and broadcasts, so every rank must call it.
double read_start_time_from_wrfinput (int lev, const std::string& fname);
#endif

// ---------------------------------------------------------------------------
//  Calendar helpers shared with ERF_NOAHMP_Advance.cpp (declared there; they
//  belong in ERF_NOAHMP.H once that header is in scope for this change).
//
//  Noah-MP's JULIAN is 1-based day-of-year plus the UTC fraction of the day:
//  1.0 == Jan 1 00Z. That is CAL_MON_DAY's convention
//  (Submodules/Noah-MP/utility/CalMonDayMod.F90 maps JULDAY = 1 -> Jan 1) and
//  the same convention ERF_Radiation.cpp:1096 uses for `calday`, so a JULIAN
//  can be handed to orbital_decl() unchanged.
//
//  Both directions go through timegm/gmtime rather than hand-rolled leap
//  arithmetic, so multi-year runs and leap years cost nothing extra.
// ---------------------------------------------------------------------------
namespace erf_noahmp {

double
epoch_from_calendar (int yr, double julian)
{
    std::tm t0{};
    t0.tm_year = yr - 1900;
    t0.tm_mon  = 0;
    t0.tm_mday = 1;          // Jan 1 00Z of that year == JULIAN 1.0
    return static_cast<double>(timegm(&t0)) + (julian - 1.0) * 86400.0;
}

void
calendar_from_epoch (double epoch, int& yr, double& julian)
{
    const double floor_sec = std::floor(epoch);
    auto tt = static_cast<std::time_t>(floor_sec);
    std::tm ti{};
#ifdef _WIN32
    gmtime_s(&ti, &tt);
#else
    gmtime_r(&tt, &ti);
#endif
    yr     = ti.tm_year + 1900;
    julian = static_cast<double>(ti.tm_yday) + 1.0     // tm_yday is 0-based
           + (ti.tm_hour * 3600.0 + ti.tm_min * 60.0 + ti.tm_sec
              + (epoch - floor_sec)) / 86400.0;
}

} // namespace erf_noahmp

void
NOAHMP::Init (const int& lev,
              const MultiFab& cons_in,
              const Geometry& geom,
              const Geometry& geom0,
              Vector<BCRec>& domain_bcs_type,
              IntVect& refRatio,
              const Real& dt,
              Vector<Vector<std::string>>& nc_init_file)
{
    // Install Noah-MP's fatal-error handler once: route NoahmpIO_fatal() through
    // amrex::Abort so a fatal error propagates via MPI_Abort. See NoahmpFatal.H.
    static const bool noahmp_fatal_installed = []() {
        NoahmpIO_set_fatal_handler([](const char* msg){
            amrex::Abort(msg ? msg : "Noah-MP fatal error");
        });
        return true;
    }();
    amrex::ignore_unused(noahmp_fatal_installed);

    m_lev   = lev;
    m_dt    = dt;
    m_geom  = geom;
    m_geom0 = geom0;
    m_domain_bcs_type = domain_bcs_type;
    m_refRatio = refRatio;

    Box domain = geom.Domain();

    // Resolve NSOIL from erf.lsm_nsoil before building the collective LSM fabs (same
    // value the parent used via Lsm_Data_Size()); namelist NSOIL asserted below.
    m_ensure_nsoil_resolved();

    // The fixed 2D fields are identity-mapped; their names mirror the enum order.
    LsmDataMap.resize(m_lsm_data_size);
    LsmDataName.resize(m_lsm_data_size);
    for (int i(0); i < LsmData_NOAHMP::NumVars; ++i) { LsmDataMap[i] = i; }
    {
        // Names from the same registry as the enum, so they cannot drift.
        const std::vector<std::string> fixed_names = {
            NOAHMP_LSMDATA_FIELDS(NOAHMP_QUOTE)
        };
        AMREX_ALWAYS_ASSERT(int(fixed_names.size()) == LsmData_NOAHMP::NumVars);
        for (int i(0); i < LsmData_NOAHMP::NumVars; ++i) { LsmDataName[i] = fixed_names[i]; }
    }
    // Per-layer soil profile: 3 groups of m_nsoil, layer index 1-based (WRF SMOIS_k).
    {
        const char* group[m_num_soil_groups] = {"smois", "sh2o", "tslb"};
        for (int g(0); g < m_num_soil_groups; ++g) {
            for (int k(0); k < m_nsoil; ++k) {
                int idx = soil_data_idx(g,k);
                LsmDataMap[idx]  = idx;
                LsmDataName[idx] = std::string(group[g]) + "_" + std::to_string(k+1);
            }
        }
    }

    LsmFluxMap  = {LsmFlux_NOAHMP::t_flux         , LsmFlux_NOAHMP::q_flux         ,
                  LsmFlux_NOAHMP::tau13          , LsmFlux_NOAHMP::tau23          };
    LsmFluxName = {"t_flux"         , "q_flux"         ,
                   "tau13"          , "tau23"          };

    // NOTE: relies on all boxes in ba spanning zlo..zhi; otherwise dm/ba no longer
    //       line up and lsm data/flux vars can't be copied directly in a parfor.

    // Set 2D box array for lsm data
    IntVect ng(1,1,0);
    BoxArray ba = cons_in.boxArray();
    DistributionMapping dm = cons_in.DistributionMap();
    BoxList bl_lsm = ba.boxList();
    for (auto& b : bl_lsm) { b.setRange(2,0); }
    BoxArray ba_lsm(std::move(bl_lsm));

    // Set up lsm geometry
    const RealBox& dom_rb = m_geom.ProbDomain();
    const Real*    dom_dx = m_geom.CellSize();
    RealBox lsm_rb = dom_rb;
    Real lsm_dx[AMREX_SPACEDIM] = {AMREX_D_DECL(dom_dx[0],dom_dx[1],m_dz_lsm)};
    Real lsm_z_hi = dom_rb.lo(2);
    Real lsm_z_lo = lsm_z_hi - Real(m_nz_lsm)*lsm_dx[2];
    lsm_rb.setHi(2,lsm_z_hi); lsm_rb.setLo(2,lsm_z_lo);
    m_lsm_geom.define( ba_lsm.minimalBox(), lsm_rb, m_geom.Coord(), m_geom.isPeriodic() );

    // Create the data (CC), runtime-sized (fixed 2D fields + 3*m_nsoil) so soil scales
    // with NSOIL. lsm_lev0_data pointers are populated later by the parent.
    lsm_fab_data.resize(m_lsm_data_size);
    lsm_lev0_data.resize(m_lsm_data_size, nullptr);
    for (auto ivar = 0; ivar < m_lsm_data_size; ++ivar) {
        lsm_fab_data[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, ng);
        lsm_fab_data[ivar]->setVal(lsm_undefined);
    }

    // Create the fluxes (CC with ghost cells for averaging)
    for (auto ivar = 0; ivar < LsmFlux_NOAHMP::NumVars; ++ivar) {
        lsm_fab_flux[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, ng);
        lsm_fab_flux[ivar]->setVal(lsm_undefined);
    }

    m_has_nc_file = (!nc_init_file[lev].empty());
    if (m_has_nc_file) {
        Print() << "Noah-MP initialization started" << std::endl;

        // Size noahmpio_vect to the local boxes. A rank owning no boxes leaves it
        // empty and relies on the class-level m_itimestep/m_dtbl instead.
        if (cons_in.local_size() > 0) {
            noahmpio_vect.resize(cons_in.local_size(), lev);
        }

        // Pinned buffer space for all the boxes
        noahmp_input_tmp.resize(cons_in.local_size());
        noahmp_output_tmp.resize(cons_in.local_size());

        int klo = domain.smallEnd(2);

        // -------------------------------------------------------------------
        // Calendar origin (plan-slucm-coupling.md section 6, Phase 0).
        //
        // Noah-MP's YR/JULIAN used to be hard-coded to 2000-01-01 and never
        // advanced. They are coupled scalars now, and ERF -- which owns the
        // clock the atmosphere and RRTMGP already run on -- is authoritative:
        // resolve ERF's simulation start once here, pin it into each block as
        // the origin, and let Advance_With_State march JULIAN from it.
        //
        // Resolution mirrors ERF.cpp:2344-2400 exactly: the inputs-file
        // `start_datetime` if the user gave one, otherwise SIMULATION START DATE
        // from the level-0 wrfinput (ERF aborts if the two disagree, so either
        // is the same number). Level 0 is deliberate -- ERF::start_time, which
        // `elapsed_time` below is measured from, is the level-0 file's date.
        //
        // Every rank takes the same branch (same ParmParse, same nc_init_file),
        // so the collective read inside read_start_time_from_wrfinput is safe.
        // -------------------------------------------------------------------
        bool   have_erf_epoch  = false;
        double erf_start_epoch = 0.0;
        {
            const std::string datetime_format = "%Y-%m-%d %H:%M:%S";
            ParmParse pp_no_prefix("");
            std::string start_datetime;
            if (pp_no_prefix.query("start_datetime", start_datetime)) {
                if (start_datetime.length() == 16) { start_datetime += ":00"; }
                const auto t = getEpochTime(start_datetime, datetime_format);
                if (t != static_cast<std::time_t>(-1)) {
                    erf_start_epoch = static_cast<double>(t);
                    have_erf_epoch  = true;
                }
            }
#ifdef ERF_USE_NETCDF
            else if (!nc_init_file.empty() && !nc_init_file[0].empty()) {
                erf_start_epoch = read_start_time_from_wrfinput(0, nc_init_file[0][0]);
                have_erf_epoch  = true;
            }
#endif
            if (!have_erf_epoch) {
                Print() << "Noah-MP: ERF has no simulation start date; falling back to "
                           "namelist.erf START_YEAR/MONTH/DAY for the calendar origin."
                        << std::endl;
            }
        }

        // Iterate over the multifab and noahmpio objects together, using the
        // multifab to set the per-box bounds on each noahmpio object.
        int idb = 0;
        for (MFIter mfi(cons_in); mfi.isValid(); ++mfi, ++idb) {

            Box bx = mfi.tilebox();

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bx.smallEnd(2) == klo,
                "NoahMP Init: box does not start at klo; z-decomposed grids are unsupported.");

            bx.makeSlab(2,klo);

            // Pinned buffers per box; output carries the 2D outputs + 3 soil groups.
            noahmp_input_tmp[idb]  = std::make_unique<FArrayBox>(bx, NoahmpInputComp::NumComps , The_Pinned_Arena());
            noahmp_output_tmp[idb] = std::make_unique<FArrayBox>(bx, NoahmpOutputComp::NumComps + m_num_soil_groups*m_nsoil, The_Pinned_Arena());

            NoahmpIO_type* noahmpio = &noahmpio_vect[idb];

            noahmpio->blkid = idb;
            noahmpio->level = lev;
            noahmpio->ScalarInitDefault();
            noahmpio->rank = ParallelDescriptor::MyProc();
            noahmpio->comm = MPI_Comm_c2f(ParallelDescriptor::Communicator());

            // namelist.erf holds noahmpio-specific parameters, read Fortran-side.
            noahmpio->ReadNamelist();

            // Pin the calendar. NoahmpReadNamelist has just filled
            // START_YR/START_JULIAN from namelist.erf's START_YEAR/MONTH/DAY/
            // HOUR/MIN, but only because they still held the -9999 sentinel;
            // ERF's clock overrides it. YR/JULIAN are seeded to the origin so
            // InitMain() (phenology, urban LAI branching) sees a real date on
            // the very first call, before Advance_With_State has run.
            // DECLIN is deliberately left at its sentinel: only the SLUCM tile
            // reads it, and it aborts on an unstaged calendar rather than
            // casting shadows for a made-up date.
            if (have_erf_epoch) {
                int    erf_yr     = 0;
                double erf_julian = 0.0;
                erf_noahmp::calendar_from_epoch(erf_start_epoch, erf_yr, erf_julian);

                if (idb == 0) {
                    const double nml_epoch = erf_noahmp::epoch_from_calendar(
                        noahmpio->START_YR, static_cast<double>(noahmpio->START_JULIAN));
                    if (std::abs(nml_epoch - erf_start_epoch) > 1.0) {
                        Print() << "***** WARNING: namelist.erf start date "
                                << getTimestamp(nml_epoch, "%Y-%m-%d %H:%M:%S", false)
                                << " UTC does not match ERF's simulation start "
                                << getTimestamp(erf_start_epoch, "%Y-%m-%d %H:%M:%S", false)
                                << " UTC. ERF's clock wins, so Noah-MP phenology and the "
                                   "SLUCM solar geometry follow ERF; fix START_YEAR/"
                                   "START_MONTH/START_DAY/START_HOUR/START_MIN in "
                                   "namelist.erf." << std::endl;
                    }
                }

                noahmpio->START_YR     = erf_yr;
                noahmpio->START_JULIAN = static_cast<noahmp_real>(erf_julian);
                noahmpio->YR           = erf_yr;
                noahmpio->JULIAN       = static_cast<noahmp_real>(erf_julian);
            }

            // Whichever branch filled it, the origin must be a real date --
            // Advance_With_State reconstructs the origin epoch from it every
            // firing, and a leftover -9999 sentinel would silently put the land
            // model tens of centuries away from the atmosphere.
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                noahmpio->START_JULIAN > noahmp_real(0.0) && noahmpio->START_YR > 0,
                "Noah-MP calendar origin was never set; namelist.erf needs "
                "START_YEAR/START_MONTH/START_DAY");

            // Assert namelist NSOIL matches erf.lsm_nsoil (used to size the fabs) so a
            // mismatch fails loudly rather than truncating the soil diagnostics.
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(noahmpio->nsoil == m_nsoil,
                "namelist.erf NSOIL does not match erf.lsm_nsoil (default 4); "
                "set erf.lsm_nsoil to the Noah-MP soil-layer count");

            // NetCDF land-file headers (also Fortran-side)
            noahmpio->ReadLandHeader();

            // Set domain/memory/tile bounds from the tile. All three are set to the
            // same bounds for now; may change for special memory management later.
            noahmpio->xstart = bx.smallEnd(0);
            noahmpio->xend   = bx.bigEnd(0);
            noahmpio->ystart = bx.smallEnd(1);
            noahmpio->yend   = bx.bigEnd(1);

            // Domain, tile, and memory bounds are all equal for now (single-slab model).
            auto set_grid_bounds = [](NoahmpIO_type* io, int x0, int x1, int y0, int y1) {
                io->ids=io->its=io->ims=x0; io->ide=io->ite=io->ime=x1;
                io->jds=io->jts=io->jms=y0; io->jde=io->jte=io->jme=y1;
                io->kds=io->kts=io->kms=1;  io->kde=io->kte=io->kme=2;
            };
            set_grid_bounds(noahmpio, noahmpio->xstart, noahmpio->xend, noahmpio->ystart, noahmpio->yend);

            // Allocate Fortran IO memory from the bounds above + namelist/header info
            noahmpio->VarInitDefault();

            // NoahmpTable.TBL input
            noahmpio->ReadTable();

            // Read/initialize from the NetCDF land file
            noahmpio->ReadLandMain();

            // Compute initial values not supplied by the land file
            noahmpio->InitMain();
        }

        // Initial land plotfile (tag 0); must run on the land-owning subcomm, not the
        // full comm, or a land-free rank never joins the collective nf90_create.
        Print() << "Noah-MP writing lnd.nc file at lev: " << lev << std::endl;
        with_land_comm([](NoahmpIO_type& noahmpio) { noahmpio.WriteLand(0); });

        // Broadcast DTBL and the initial substep counter so the firing decision is
        // identical on every rank. Land-free ranks use max-reduction-losing sentinels.
        m_dtbl      = noahmpio_vect.empty() ? std::numeric_limits<Real>::lowest()
                                            : static_cast<Real>(noahmpio_vect[0].DTBL);
        m_itimestep = noahmpio_vect.empty() ? std::numeric_limits<int>::lowest()
                                            : noahmpio_vect[0].itimestep;
        ParallelDescriptor::ReduceRealMax(m_dtbl);
        ParallelDescriptor::ReduceIntMax(m_itimestep);

        // Guard against a decomposition in which no rank owns a land box.
        AMREX_ALWAYS_ASSERT(m_dtbl > Real(0.0));
        AMREX_ALWAYS_ASSERT(m_dt <= m_dtbl);

        Print() << "Noah-MP initialization completed" << std::endl;
    } // has nc_init_file

};
