# Single-column SLUCM harness: build the SAME driver twice, once against ERF's
# vendored module_sf_urban.F90 and once against the pristine upstream WRF
# phys/module_sf_urban.F, then compare the two trajectories.
#
#   make            # build both (upstream target is skipped if absent)
#   make run        # run both and compare
#   make erf        # ERF side only -- works with no WRF source at all
#   make trap       # rebuild with FPE trapping + snan init, run both
#   make clean
#
# No NetCDF, no MPI, no AMReX: plain gfortran. ERF proper cannot be built on
# macOS because ERF_NCInterface.H needs a parallel netcdf-c, but nothing in this
# harness includes it, so a laptop is a perfectly good host.

FC      = gfortran

NOAHMP  := ../../../Submodules/Noah-MP
ERFDRV  := $(NOAHMP)/drivers/erf

# The vendored file is bare REAL and ERF compiles it with real promotion; the
# driver asserts at run time that this actually took effect. Both sides get the
# identical flag set -- that is what makes the comparison meaningful.
PROMOTE := -fdefault-real-8 -fdefault-double-8
FFLAGS  = -O2 -g -ffree-form -ffree-line-length-none -cpp $(PROMOTE)

# `make trap` adds these. The 3-D smoke test cannot trap FPEs at all -- Noah-MP
# writes through HDF5, and H5T__init_native_float_types raises FE_INVALID while
# probing native float types during library init. Offline there is no HDF5, so
# this harness can trap, which is how the upstream undefined reads become
# visible rather than merely argued.
TRAPFLAGS := -ffpe-trap=invalid,zero,overflow -finit-real=snan -fbacktrace

UPSTREAM ?= upstream/module_sf_urban.F

BUILD_ERF := build-erf
BUILD_WRF := build-wrf

.PHONY: all erf wrf run trap clean

all: erf wrf

# ---------------------------------------------------------------- ERF side ---
erf: slucm_column_erf

slucm_column_erf: $(ERFDRV)/module_sf_urban.F90 $(ERFDRV)/NoahmpUrbanShimMod.F90 \
                  $(NOAHMP)/utility/Machine.F90 slucm_column.F90
	@mkdir -p $(BUILD_ERF)
	$(FC) $(FFLAGS) -J$(BUILD_ERF) -c $(NOAHMP)/utility/Machine.F90     -o $(BUILD_ERF)/Machine.o
	$(FC) $(FFLAGS) -J$(BUILD_ERF) -c $(ERFDRV)/NoahmpUrbanShimMod.F90  -o $(BUILD_ERF)/shim.o
	$(FC) $(FFLAGS) -J$(BUILD_ERF) -c $(ERFDRV)/module_sf_urban.F90     -o $(BUILD_ERF)/urban.o
	$(FC) $(FFLAGS) -DHAVE_KIND_URBAN_REAL -I$(BUILD_ERF) -J$(BUILD_ERF) -c slucm_column.F90   -o $(BUILD_ERF)/driver.o
	$(FC) $(FFLAGS) -o $@ $(BUILD_ERF)/Machine.o $(BUILD_ERF)/shim.o \
	                      $(BUILD_ERF)/urban.o $(BUILD_ERF)/driver.o

# ----------------------------------------------------------- upstream side ---
# The upstream file is NEVER edited. It is .F (fixed-form by extension) but
# free-form in content, so -ffree-form is mandatory; that is exactly why the
# vendored copy was renamed to .F90.
wrf:
	@if [ -f "$(UPSTREAM)" ]; then \
	   $(MAKE) slucm_column_wrf; \
	 else \
	   echo "note: $(UPSTREAM) not present -- skipping the upstream build."; \
	   echo "      fetch it with: make fetch-upstream"; \
	 fi

slucm_column_wrf: $(UPSTREAM) wrf_stubs.F90 slucm_column.F90
	@mkdir -p $(BUILD_WRF)
	$(FC) $(FFLAGS) -J$(BUILD_WRF) -c wrf_stubs.F90      -o $(BUILD_WRF)/stubs.o
	$(FC) $(FFLAGS) -J$(BUILD_WRF) -c $(UPSTREAM)        -o $(BUILD_WRF)/urban.o
	$(FC) $(FFLAGS) -I$(BUILD_WRF) -J$(BUILD_WRF) -c slucm_column.F90 -o $(BUILD_WRF)/driver.o
	$(FC) $(FFLAGS) -o $@ $(BUILD_WRF)/stubs.o $(BUILD_WRF)/urban.o $(BUILD_WRF)/driver.o

# Recover the exact revision the vendored copy came from, identified by md5 in
# the vendored file's provenance header.
.PHONY: fetch-upstream
fetch-upstream:
	@mkdir -p upstream
	@echo "fetching wrf-model/WRF @ 8fa379b4 : phys/module_sf_urban.F"
	@rm -rf upstream/.wrf
	@git clone --filter=blob:none --no-checkout --quiet \
	     https://github.com/wrf-model/WRF upstream/.wrf
	@git -C upstream/.wrf show \
	     8fa379b458415a120fa8328cff78194b42c64a35:phys/module_sf_urban.F \
	     > $(UPSTREAM)
	@rm -rf upstream/.wrf
	@printf 'expected md5 ad16e13c8f24b30c9fa815fe10702db1\n     got md5 '
	@md5 -q $(UPSTREAM) 2>/dev/null || md5sum $(UPSTREAM) | cut -d' ' -f1

# ------------------------------------------------------------------- run -----
run: all URBPARM.TBL
	./slucm_column_erf forcing_48h.csv trajectory_erf.csv
	@if [ -x ./slucm_column_wrf ]; then \
	   ./slucm_column_wrf forcing_48h.csv trajectory_wrf.csv; \
	   ./compare_trajectories.py trajectory_erf.csv trajectory_wrf.csv; \
	 else \
	   echo "upstream build absent -- ERF trajectory only, no comparison made."; \
	 fi

trap:
	$(MAKE) clean
	$(MAKE) FFLAGS="-O0 -g -ffree-form -ffree-line-length-none -cpp $(PROMOTE) $(TRAPFLAGS)" run

URBPARM.TBL:
	ln -sf $(ERFDRV)/URBPARM.TBL URBPARM.TBL

clean:
	rm -rf $(BUILD_ERF) $(BUILD_WRF) slucm_column_erf slucm_column_wrf \
	       trajectory_erf.csv trajectory_wrf.csv
