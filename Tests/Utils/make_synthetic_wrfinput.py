#!/usr/bin/env python3
"""Generate a tiny synthetic ``wrfinput`` NetCDF file that exercises ERF's
real-data (``erf.init_type = "WRFInput"``) path *including* the urban tile.

Why this exists
---------------
``Exec/RegTests/WPS_Test`` needs ``wrfinput_chisholmview_d01`` (and a matching
``wrfbdy``), neither of which is in the repository, and the ERF gold-file
repository is not reachable from GitHub Actions.  Phase 1 item 2 of
``Source/LandSurfaceModel/Noah-MP/dev/plan-slucm-coupling.md`` is explicitly
blocked on that data decision.  This script sidesteps it for the *urban*
purpose: it fabricates a 10x10 domain small enough to commit or to regenerate
in CI, containing every field that
``Source/Initialization/ERF_InitFromWRFInput.cpp``,
``Source/IO/ERF_ReadFromWRFInput.cpp`` and
``Submodules/Noah-MP/drivers/erf/NoahmpReadLandMod.F90`` read, plus the urban
statics (``FRC_URB2D``, ``IVGTYP == ISURBAN``) that the SLUCM tile keys off.

No boundary file is produced.  Run with ``erf.real_width = 0`` and no
``erf.nc_bdy_file`` so nothing needs a fabricated ``wrfbdy``.

What the readers require (and this script honours)
-------------------------------------------------
* Every variable ERF reads must have ``Time`` as its first dimension --
  ``ERF_NCWpsFile.H`` asserts ``dimnames[0] == "Time" || "time"``.
* The dimensions ``west_east_stag``, ``south_north_stag`` and ``bottom_top``
  must exist; ``read_from_wrfinput`` queries them by name and asserts
  ``west_east_stag == WEST-EAST_GRID_DIMENSION``.
* Density: ``CheckForDensity`` (``ERF_InitFromWRFInput.cpp``) accepts *either*
  ``ALB`` + ``AL`` *or* ``ALT`` alone, and aborts when both are present.  A real
  WRF ``wrfinput`` carries all three, so the default here is ``ALB`` + ``AL``;
  use ``--density alt`` for the other branch.
* ``IVGTYP``/``ISLTYP``/``UTYPE_URB2D`` are written as ``NC_INT`` (as WRF does).
  Noah-MP reads them into a ``real`` array and ERF into a ``float`` buffer;
  NetCDF converts in both cases.
* ``MMINLU`` selects the USGS vs MODIS block of ``NoahmpTable.TBL``, and the
  ``ISURBAN`` global attribute must match that block's ``ISURBAN``
  (USGS 1, MODIS 13) or no cell is recognised as urban.

Fields that are deliberately *not* written: ``LF_URB2D`` (4-D, per-canyon-
direction) and ``URB_PARAM`` (the WPS/NUDAPT stack).  ``urban_var_init``
falls back to ``URBPARM.TBL`` geometry whenever ``HGT_URB2D <= 0``, so Phase 1
does not need them; pass ``--no-urban-morphology`` to drop the optional 2-D
morphology fields as well and rely purely on the table fallback.

Usage
-----
    python3 Tests/Utils/make_synthetic_wrfinput.py -o wrfinput_urban_d01
    python3 Tests/Utils/make_synthetic_wrfinput.py -o out.nc --dump

Requires ``netCDF4`` (which brings ``numpy``).  Nothing else.
"""

from __future__ import annotations

import argparse
import math
import sys

try:
    import numpy as np
except ImportError:  # pragma: no cover - trivial guard
    sys.exit("ERROR: numpy is required (it ships with netCDF4): pip install netCDF4")

try:
    from netCDF4 import Dataset
except ImportError:  # pragma: no cover - trivial guard
    sys.exit("ERROR: the netCDF4 module is required: pip install netCDF4")


# --------------------------------------------------------------------------
# Physical constants.  These only need to be mutually consistent *within this
# file* -- the point is a hydrostatic, invertible base state, not agreement
# with ERF_Constants.H to the last digit.
# --------------------------------------------------------------------------
G = 9.81  # m s-2
R_D = 287.0  # J kg-1 K-1
CP = 1004.5  # J kg-1 K-1
P_0 = 1.0e5  # Pa
T_SFC = 290.0  # K, temperature at z = 0
LAPSE = 0.0065  # K m-1
THETA_REF = 300.0  # K, the constant ERF adds back to THM

# Land-use table indices, from Submodules/Noah-MP/drivers/erf/NoahmpTable.TBL.
LANDUSE_TABLES = {
    # name: (ISURBAN, ISWATER, ISICE, NATURAL, NUM_LAND_CAT)
    "USGS": (1, 16, 24, 5, 24),
    "MODIFIED_IGBP_MODIS_NOAH": (13, 17, 15, 14, 20),
}

def _cdl_type(dtype) -> str:
    return {"f4": "float", "f8": "double", "i4": "int", "S1": "char"}.get(
        np.dtype(dtype).str[1:], str(dtype)
    )


class WrfInputBuilder:
    """Thin wrapper that records everything written so ``--dump`` can replay it."""

    def __init__(self, path: str, fmt: str):
        self.nc = Dataset(path, "w", format=fmt)
        self.path = path

    def dim(self, name: str, size):
        self.nc.createDimension(name, size)

    def gattr(self, name: str, value):
        # netCDF4 maps a Python int to NC_INT64 unless we hand it a numpy
        # scalar; WRF (and ERF's nc_get_att_int) expect NC_INT / NC_FLOAT.
        self.nc.setncattr(name, value)

    @staticmethod
    def _memory_order(dims) -> str:
        spatial = [d for d in dims if d != "Time"]
        has_x = any(d.startswith("west_east") for d in spatial)
        has_y = any(d.startswith("south_north") for d in spatial)
        has_z = any(d.startswith(("bottom_top", "soil_layers")) for d in spatial)
        if has_x and has_y and has_z:
            return "XYZ"
        if has_x and has_y:
            return "XY "
        if has_z:
            return "Z  "
        return "0  "

    @staticmethod
    def _stagger(dims) -> str:
        s = ""
        if any(d == "west_east_stag" for d in dims):
            s += "X"
        if any(d == "south_north_stag" for d in dims):
            s += "Y"
        if any(d in ("bottom_top_stag", "soil_layers_stag") for d in dims):
            s += "Z"
        return s

    def var(self, name, dtype, dims, data, description="", units=""):
        v = self.nc.createVariable(name, dtype, dims)
        v.FieldType = {"f4": 104, "i4": 106, "S1": 106}[np.dtype(dtype).str[1:]]
        v.MemoryOrder = self._memory_order(dims)
        v.description = description
        v.units = units
        v.stagger = self._stagger(dims)
        v[...] = data
        return v

    def close(self):
        self.nc.close()


def build(args) -> str:
    nx, ny, nz, nsoil = args.nx, args.ny, args.nz, args.nsoil
    dx, dy, dz = args.dx, args.dy, args.dz

    if args.mminlu not in LANDUSE_TABLES:
        sys.exit(f"ERROR: unknown MMINLU '{args.mminlu}'; "
                 f"known: {sorted(LANDUSE_TABLES)}")
    isurban, iswater, isice, natural, num_land_cat = LANDUSE_TABLES[args.mminlu]

    # ---------------- vertical / thermodynamic base state ----------------
    # Staggered (w) heights and unstaggered (mass) heights.
    z_w = dz * np.arange(nz + 1, dtype=np.float64)
    z_c = 0.5 * (z_w[:-1] + z_w[1:])

    def p_of_z(z):
        """Hydrostatic pressure for a constant-lapse-rate atmosphere."""
        return P_0 * (1.0 - LAPSE * z / T_SFC) ** (G / (R_D * LAPSE))

    t_c = T_SFC - LAPSE * z_c  # K
    pb_c = p_of_z(z_c)  # base pressure on mass levels, Pa
    pb_w = p_of_z(z_w)
    p_top = float(pb_w[-1])

    theta = t_c * (P_0 / pb_c) ** (R_D / CP)
    thm = theta - THETA_REF  # WRF stores the perturbation from 300 K
    alb = R_D * t_c / pb_c  # 1/rho for the base state
    al = np.zeros_like(alb)  # perturbation density is zero: P == 0
    alt = alb + al

    phb = G * z_w  # base geopotential
    ph = np.zeros_like(phb)

    # WRF mass coordinate.  eta = 1 at the surface, 0 at the model top.
    mub = P_0 - p_top
    eta_w = (pb_w - p_top) / mub
    dnw = eta_w[1:] - eta_w[:-1]  # negative, k increasing upward
    rdnw = 1.0 / dnw
    # Pure sigma (non-hybrid) reduction of WRF's hybrid coordinate.
    c1h = np.ones(nz)
    c2h = np.zeros(nz)

    # ---------------- horizontal / surface fields ----------------
    lat0, lon0 = args.lat, args.lon
    dlat = dy / 111000.0
    dlon = dx / (111000.0 * math.cos(math.radians(lat0)))

    j_c = (np.arange(ny) - 0.5 * (ny - 1))
    i_c = (np.arange(nx) - 0.5 * (nx - 1))
    j_s = (np.arange(ny + 1) - 0.5 * ny)
    i_s = (np.arange(nx + 1) - 0.5 * nx)

    xlat = lat0 + dlat * j_c[:, None] * np.ones((1, nx))
    xlong = lon0 + dlon * np.ones((ny, 1)) * i_c[None, :]
    xlat_v = lat0 + dlat * j_s[:, None] * np.ones((1, nx))  # (south_north_stag, west_east)
    xlong_u = lon0 + dlon * np.ones((ny, 1)) * i_s[None, :]  # (south_north, west_east_stag)

    # ---------------- the urban block ----------------
    i0, i1 = args.urban_i0, args.urban_i1
    j0, j1 = args.urban_j0, args.urban_j1
    if not (0 <= i0 <= i1 < nx and 0 <= j0 <= j1 < ny):
        sys.exit("ERROR: --urban-i0/--urban-i1/--urban-j0/--urban-j1 out of range")

    urban = np.zeros((ny, nx), dtype=bool)
    urban[j0:j1 + 1, i0:i1 + 1] = True
    n_urban = int(urban.sum())
    if n_urban == 0:
        sys.exit("ERROR: refusing to write a file with zero urban cells")

    ivgtyp = np.full((ny, nx), args.rural_vegtype, dtype=np.int32)
    ivgtyp[urban] = isurban
    isltyp = np.full((ny, nx), args.soiltype, dtype=np.int32)

    frc_urb = np.zeros((ny, nx), dtype=np.float64)
    frc_urb[urban] = args.urban_fraction

    # urban_var_init assigns UTYPE_URB2D = 2 ("high-density residential") to
    # every ISURBAN cell when use_wudapt_lcz = 0; write the same thing so the
    # file is self-describing rather than relying on that fallback.
    utype_urb = np.zeros((ny, nx), dtype=np.int32)
    utype_urb[urban] = args.urban_type

    def urban_field(urban_value, rural_value):
        a = np.full((ny, nx), rural_value, dtype=np.float64)
        a[urban] = urban_value
        return a

    tsk = urban_field(args.tsk_urban, args.tsk_rural)
    vegfra = urban_field(10.0, 50.0)
    lai = urban_field(0.5, 2.0)
    shdmin = urban_field(5.0, 5.0)
    # Deliberately low over the city: Noah-MP's rural-tile contract
    # (ConfigVarInTransferMod.F90) overrides GVFMAX to 96 % when
    # SF_URBAN_PHYSICS > 0, and that override is only observable if the file
    # disagrees with it.
    shdmax = urban_field(10.0, 80.0)

    ones2d = np.ones((ny, nx))
    zeros2d = np.zeros((ny, nx))

    # Soil: WRF's default 4-layer configuration.
    if nsoil == 4:
        dzs = np.array([0.10, 0.30, 0.60, 1.00])
        zs = np.array([0.05, 0.25, 0.70, 1.50])
    else:
        dzs = np.full(nsoil, 2.0 / nsoil)
        zs = np.cumsum(dzs) - 0.5 * dzs
    tslb = (290.0 - np.arange(nsoil))[:, None, None] * np.ones((1, ny, nx))
    smois = (0.25 + 0.01 * np.arange(nsoil))[:, None, None] * np.ones((1, ny, nx))

    def col3d(profile, shape_yx=(ny, nx)):
        return profile[:, None, None] * np.ones((1,) + shape_yx)

    # ---------------- write ----------------
    b = WrfInputBuilder(args.output, args.format)

    b.dim("Time", None)
    b.dim("DateStrLen", 19)
    b.dim("west_east", nx)
    b.dim("south_north", ny)
    b.dim("west_east_stag", nx + 1)
    b.dim("south_north_stag", ny + 1)
    b.dim("bottom_top", nz)
    b.dim("bottom_top_stag", nz + 1)
    b.dim("soil_layers_stag", nsoil)

    i4 = lambda v: np.int32(v)      # noqa: E731 - keep NC_INT, not NC_INT64
    f4 = lambda v: np.float32(v)    # noqa: E731 - keep NC_FLOAT, not NC_DOUBLE

    b.gattr("TITLE", " OUTPUT FROM ERF Tests/Utils/make_synthetic_wrfinput.py")
    b.gattr("SIMULATION_START_DATE", args.start_date)
    b.gattr("WEST-EAST_GRID_DIMENSION", i4(nx + 1))
    b.gattr("SOUTH-NORTH_GRID_DIMENSION", i4(ny + 1))
    b.gattr("BOTTOM-TOP_GRID_DIMENSION", i4(nz + 1))
    b.gattr("DX", f4(dx))
    b.gattr("DY", f4(dy))
    b.gattr("DT", f4(args.dt))
    b.gattr("GRID_ID", i4(1))
    b.gattr("PARENT_ID", i4(0))
    b.gattr("I_PARENT_START", i4(1))
    b.gattr("J_PARENT_START", i4(1))
    b.gattr("PARENT_GRID_RATIO", i4(1))
    b.gattr("CEN_LAT", f4(lat0))
    b.gattr("CEN_LON", f4(lon0))
    b.gattr("TRUELAT1", f4(args.truelat1))
    b.gattr("TRUELAT2", f4(args.truelat2))
    b.gattr("MOAD_CEN_LAT", f4(lat0))
    b.gattr("STAND_LON", f4(lon0))
    b.gattr("MAP_PROJ", i4(1))
    b.gattr("MAP_PROJ_CHAR", "Lambert Conformal")
    b.gattr("MMINLU", args.mminlu)
    b.gattr("NUM_LAND_CAT", i4(num_land_cat))
    b.gattr("ISWATER", i4(iswater))
    b.gattr("ISLAKE", i4(-1))
    b.gattr("ISICE", i4(isice))
    b.gattr("ISURBAN", i4(isurban))
    b.gattr("ISOILWATER", i4(14))
    b.gattr("USE_THETA_M", i4(0))
    b.gattr("P_TOP", f4(p_top))
    # The plan's Phase 0 wants JULIAN/YR/GMT wired through the ABI; carry them
    # here so the file is not the thing standing in the way.
    b.gattr("JULYR", i4(args.julyr))
    b.gattr("JULDAY", i4(args.julday))
    b.gattr("GMT", f4(args.gmt))

    times = np.array([list(args.start_date.ljust(19)[:19].encode("ascii"))],
                     dtype="S1").reshape(1, 19)
    tv = b.nc.createVariable("Times", "S1", ("Time", "DateStrLen"))
    tv[...] = times

    SN_WE = ("Time", "south_north", "west_east")
    SNs_WE = ("Time", "south_north_stag", "west_east")
    SN_WEs = ("Time", "south_north", "west_east_stag")
    BT_SN_WE = ("Time", "bottom_top", "south_north", "west_east")
    BTs_SN_WE = ("Time", "bottom_top_stag", "south_north", "west_east")
    BT_SN_WEs = ("Time", "bottom_top", "south_north", "west_east_stag")
    BT_SNs_WE = ("Time", "bottom_top", "south_north_stag", "west_east")
    SL_SN_WE = ("Time", "soil_layers_stag", "south_north", "west_east")

    # --- dycore 3-D state -------------------------------------------------
    b.var("U", "f4", BT_SN_WEs,
          np.full((1, nz, ny, nx + 1), args.u), "x-wind component", "m s-1")
    b.var("V", "f4", BT_SNs_WE,
          np.full((1, nz, ny + 1, nx), args.v), "y-wind component", "m s-1")
    b.var("W", "f4", BTs_SN_WE,
          np.zeros((1, nz + 1, ny, nx)), "z-wind component", "m s-1")
    b.var("THM", "f4", BT_SN_WE, col3d(thm)[None, ...],
          "perturbation potential temperature (theta - t0)", "K")
    b.var("PH", "f4", BTs_SN_WE, col3d(ph)[None, ...],
          "perturbation geopotential", "m2 s-2")
    b.var("PHB", "f4", BTs_SN_WE, col3d(phb)[None, ...],
          "base-state geopotential", "m2 s-2")
    b.var("P", "f4", BT_SN_WE, np.zeros((1, nz, ny, nx)),
          "perturbation pressure", "Pa")
    b.var("PB", "f4", BT_SN_WE, col3d(pb_c)[None, ...],
          "BASE STATE PRESSURE", "Pa")

    # CheckForDensity() aborts if ALB/AL and ALT are both present.
    if args.density == "alb_al":
        b.var("ALB", "f4", BT_SN_WE, col3d(alb)[None, ...],
              "inverse base density", "m3 kg-1")
        b.var("AL", "f4", BT_SN_WE, col3d(al)[None, ...],
              "inverse perturbation density", "m3 kg-1")
    else:
        b.var("ALT", "f4", BT_SN_WE, col3d(alt)[None, ...],
              "inverse density", "m3 kg-1")

    if args.moist:
        qv = 0.008 * np.exp(-z_c / 3000.0)
        b.var("QVAPOR", "f4", BT_SN_WE, col3d(qv)[None, ...],
              "Water vapor mixing ratio", "kg kg-1")
        b.var("QCLOUD", "f4", BT_SN_WE, np.zeros((1, nz, ny, nx)),
              "Cloud water mixing ratio", "kg kg-1")
        b.var("QRAIN", "f4", BT_SN_WE, np.zeros((1, nz, ny, nx)),
              "Rain water mixing ratio", "kg kg-1")

    # --- 1-D coordinate metadata -----------------------------------------
    b.var("C1H", "f4", ("Time", "bottom_top"), c1h[None, :],
          "half levels, c1h = d(bf)/d(eta)", "Dimensionless")
    b.var("C2H", "f4", ("Time", "bottom_top"), c2h[None, :],
          "half levels, c2h = (1-c1h)*(p0-pt)", "Pa")
    b.var("RDNW", "f4", ("Time", "bottom_top"), rdnw[None, :],
          "inverse d(eta) values between full (w) levels", "")
    b.var("ZS", "f4", ("Time", "soil_layers_stag"), zs[None, :],
          "DEPTHS OF CENTERS OF SOIL LAYERS", "m")
    b.var("DZS", "f4", ("Time", "soil_layers_stag"), dzs[None, :],
          "THICKNESSES OF SOIL LAYERS", "m")

    # --- 2-D geometry / map factors ---------------------------------------
    b.var("XLAT", "f4", SN_WE, xlat[None, ...], "LATITUDE, SOUTH IS NEGATIVE",
          "degree_north")
    b.var("XLONG", "f4", SN_WE, xlong[None, ...], "LONGITUDE, WEST IS NEGATIVE",
          "degree_east")
    b.var("XLAT_V", "f4", SNs_WE, xlat_v[None, ...],
          "LATITUDE, SOUTH IS NEGATIVE", "degree_north")
    b.var("XLONG_U", "f4", SN_WEs, xlong_u[None, ...],
          "LONGITUDE, WEST IS NEGATIVE", "degree_east")
    b.var("MAPFAC_M", "f4", SN_WE, ones2d[None, ...],
          "Map scale factor on mass grid", "")
    b.var("MAPFAC_MX", "f4", SN_WE, ones2d[None, ...],
          "Map scale factor on mass grid, x direction", "")
    b.var("MAPFAC_MY", "f4", SN_WE, ones2d[None, ...],
          "Map scale factor on mass grid, y direction", "")
    b.var("MAPFAC_U", "f4", SN_WEs, np.ones((1, ny, nx + 1)),
          "Map scale factor on u-grid", "")
    b.var("MAPFAC_V", "f4", SNs_WE, np.ones((1, ny + 1, nx)),
          "Map scale factor on v-grid", "")
    b.var("MUB", "f4", SN_WE, np.full((1, ny, nx), mub),
          "base state dry air mass in column", "Pa")
    b.var("PSFC", "f4", SN_WE, np.full((1, ny, nx), P_0), "SFC PRESSURE", "Pa")

    # --- 2-D land surface -------------------------------------------------
    b.var("HGT", "f4", SN_WE, zeros2d[None, ...], "Terrain Height", "m")
    b.var("LANDMASK", "f4", SN_WE, ones2d[None, ...],
          "LAND MASK (1 FOR LAND, 0 FOR WATER)", "")
    b.var("XLAND", "f4", SN_WE, ones2d[None, ...],
          "LAND MASK (1 FOR LAND, 2 FOR WATER)", "")
    b.var("SEAICE", "f4", SN_WE, zeros2d[None, ...], "SEA ICE FLAG", "")
    b.var("TSK", "f4", SN_WE, tsk[None, ...], "SURFACE SKIN TEMPERATURE", "K")
    b.var("SST", "f4", SN_WE, np.full((1, ny, nx), 300.0),
          "SEA SURFACE TEMPERATURE", "K")
    b.var("TMN", "f4", SN_WE, np.full((1, ny, nx), 287.0),
          "SOIL TEMPERATURE AT LOWER BOUNDARY", "K")
    b.var("CANWAT", "f4", SN_WE, zeros2d[None, ...], "CANOPY WATER", "kg m-2")
    b.var("SNOW", "f4", SN_WE, zeros2d[None, ...], "SNOW WATER EQUIVALENT",
          "kg m-2")
    b.var("SNOWC", "f4", SN_WE, zeros2d[None, ...], "FLAG INDICATING SNOW COVERAGE",
          "")
    b.var("SNOWH", "f4", SN_WE, zeros2d[None, ...], "PHYSICAL SNOW DEPTH", "m")
    b.var("VEGFRA", "f4", SN_WE, vegfra[None, ...], "VEGETATION FRACTION", "%")
    b.var("LAI", "f4", SN_WE, lai[None, ...], "LEAF AREA INDEX", "m-2 m-2")
    b.var("SHDMIN", "f4", SN_WE, shdmin[None, ...],
          "ANNUAL MIN VEG FRACTION", "%")
    b.var("SHDMAX", "f4", SN_WE, shdmax[None, ...],
          "ANNUAL MAX VEG FRACTION", "%")
    b.var("IVGTYP", "i4", SN_WE, ivgtyp[None, ...], "DOMINANT VEGETATION CATEGORY",
          "")
    b.var("ISLTYP", "i4", SN_WE, isltyp[None, ...], "DOMINANT SOIL CATEGORY", "")

    # --- soil ------------------------------------------------------------
    b.var("TSLB", "f4", SL_SN_WE, tslb[None, ...], "SOIL TEMPERATURE", "K")
    b.var("SMOIS", "f4", SL_SN_WE, smois[None, ...], "SOIL MOISTURE", "m3 m-3")
    b.var("SH2O", "f4", SL_SN_WE, smois[None, ...],
          "SOIL LIQUID WATER", "m3 m-3")

    # --- urban ------------------------------------------------------------
    # FRC_URB2D is the one urban static ERF itself reads
    # (ERF_ReadFromWRFInput.cpp / ERF_InitFromWRFInput.cpp).
    b.var("FRC_URB2D", "f4", SN_WE, frc_urb[None, ...], "URBAN FRACTION", "")
    if args.urban_morphology:
        b.var("UTYPE_URB2D", "i4", SN_WE, utype_urb[None, ...], "URBAN TYPE", "")
        b.var("LP_URB2D", "f4", SN_WE, urban_field(0.35, 0.0)[None, ...],
              "PLAN AREA FRACTION", "")
        b.var("LB_URB2D", "f4", SN_WE, urban_field(1.20, 0.0)[None, ...],
              "BUILDING SURFACE AREA TO PLAN AREA RATIO", "")
        b.var("HGT_URB2D", "f4", SN_WE, urban_field(args.building_height, 0.0)[None, ...],
              "AVERAGE BUILDING HEIGHT WEIGHTED BY BUILDING PLAN AREA", "m")
        b.var("MH_URB2D", "f4", SN_WE, urban_field(args.building_height, 0.0)[None, ...],
              "MEAN BUILDING HEIGHT", "m")
        b.var("STDH_URB2D", "f4", SN_WE, urban_field(2.0, 0.0)[None, ...],
              "STANDARD DEVIATION OF BUILDING HEIGHT", "m")

    b.close()

    print(f"wrote {args.output}")
    print(f"  domain            : {nx} x {ny} x {nz}, dx = dy = {dx} m, dz = {dz} m")
    print(f"  soil layers       : {nsoil}")
    print(f"  land-use table    : {args.mminlu} (ISURBAN = {isurban}, "
          f"NATURAL = {natural})")
    print(f"  urban cells       : {n_urban} of {nx * ny} "
          f"(i {i0}..{i1}, j {j0}..{j1}), IVGTYP = {isurban}")
    print(f"  FRC_URB2D         : {args.urban_fraction} over those cells, 0 elsewhere")
    print(f"  density variables : "
          f"{'ALB + AL' if args.density == 'alb_al' else 'ALT'}")
    print(f"  p_top             : {p_top:.2f} Pa,  MUB = {mub:.2f} Pa")
    print(f"  lowest mass level : z = {z_c[0]:.2f} m "
          f"(URBPARM.TBL type 2 needs ZA > 8.05 m, and >= 2*ZR = 15 m for MOST)")
    return args.output


def dump_header(path: str) -> None:
    """An ``ncdump -h``-equivalent, so the file can be inspected without the
    NetCDF command-line tools installed."""
    nc = Dataset(path, "r")
    print(f"netcdf {path.rsplit('/', 1)[-1].rsplit('.', 1)[0]} {{")
    print("dimensions:")
    for name, d in nc.dimensions.items():
        size = "UNLIMITED ; // (%d currently)" % len(d) if d.isunlimited() else f"{len(d)} ;"
        print(f"\t{name} = {size}")
    print("variables:")
    for name, v in nc.variables.items():
        dims = ", ".join(v.dimensions)
        print(f"\t{_cdl_type(v.dtype)} {name}({dims}) ;")
        for a in v.ncattrs():
            val = v.getncattr(a)
            val = f'"{val}"' if isinstance(val, str) else val
            print(f"\t\t{name}:{a} = {val} ;")
    print("\n// global attributes:")
    for a in nc.ncattrs():
        val = nc.getncattr(a)
        val = f'"{val}"' if isinstance(val, str) else val
        print(f"\t\t:{a} = {val} ;")
    print("}")
    nc.close()


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-o", "--output", default="wrfinput_urban_d01",
                   help="output NetCDF file")
    p.add_argument("--format", default="NETCDF4",
                   choices=["NETCDF4", "NETCDF4_CLASSIC", "NETCDF3_64BIT_OFFSET",
                            "NETCDF3_CLASSIC"],
                   help="NetCDF file format")
    p.add_argument("--dump", action="store_true",
                   help="print an ncdump -h style header after writing")

    g = p.add_argument_group("grid")
    g.add_argument("--nx", type=int, default=10)
    g.add_argument("--ny", type=int, default=10)
    g.add_argument("--nz", type=int, default=20)
    g.add_argument("--nsoil", type=int, default=4,
                   help="must be 4 for SLUCM: urban_var_init hard-codes "
                        "TRL_URB3D(I,1..4,J)")
    g.add_argument("--dx", type=float, default=1000.0,
                   help="plan section 7.2 supported envelope is dx >~ 1 km")
    g.add_argument("--dy", type=float, default=1000.0)
    g.add_argument("--dz", type=float, default=100.0,
                   help="uniform layer thickness; the first mass level sits at "
                        "dz/2 and must clear 2*max(ZR)")
    g.add_argument("--lat", type=float, default=40.0)
    g.add_argument("--lon", type=float, default=-105.0)
    g.add_argument("--truelat1", type=float, default=30.0)
    g.add_argument("--truelat2", type=float, default=60.0)
    g.add_argument("--dt", type=float, default=6.0)

    t = p.add_argument_group("time")
    t.add_argument("--start-date", default="2020-07-01_12:00:00",
                   help="local solar noon-ish, so the SLUCM shadow model is "
                        "actually exercised")
    t.add_argument("--julyr", type=int, default=2020)
    t.add_argument("--julday", type=int, default=183)
    t.add_argument("--gmt", type=float, default=12.0)

    a = p.add_argument_group("atmosphere")
    a.add_argument("--u", type=float, default=5.0)
    a.add_argument("--v", type=float, default=0.0)
    a.add_argument("--density", default="alb_al", choices=["alb_al", "alt"],
                   help="which density representation to write; ERF aborts if "
                        "both are present")
    a.add_argument("--moist", action=argparse.BooleanOptionalAction, default=True,
                   help="write QVAPOR/QCLOUD/QRAIN")

    u = p.add_argument_group("land / urban")
    u.add_argument("--mminlu", default="USGS", choices=sorted(LANDUSE_TABLES))
    u.add_argument("--rural-vegtype", type=int, default=5,
                   help="USGS 5 = cropland/grassland mosaic (= NATURAL)")
    u.add_argument("--soiltype", type=int, default=8,
                   help="USGS/STATSGO 8 = silty clay loam")
    u.add_argument("--urban-fraction", type=float, default=0.9)
    u.add_argument("--urban-type", type=int, default=2,
                   help="1 low-density res., 2 high-density res., 3 commercial")
    u.add_argument("--building-height", type=float, default=7.5,
                   help="ZR for URBPARM.TBL urban type 2")
    u.add_argument("--urban-i0", type=int, default=3)
    u.add_argument("--urban-i1", type=int, default=6)
    u.add_argument("--urban-j0", type=int, default=3)
    u.add_argument("--urban-j1", type=int, default=6)
    u.add_argument("--urban-morphology",
                   action=argparse.BooleanOptionalAction, default=True,
                   help="write the optional 2-D morphology fields "
                        "(UTYPE/LP/LB/HGT/MH/STDH_URB2D); with --no-urban-"
                        "morphology urban_var_init falls back to URBPARM.TBL")
    u.add_argument("--tsk-rural", type=float, default=300.0)
    u.add_argument("--tsk-urban", type=float, default=302.0)

    args = p.parse_args(argv)
    path = build(args)
    if args.dump:
        print()
        dump_header(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
