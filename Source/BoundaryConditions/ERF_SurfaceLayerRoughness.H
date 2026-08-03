#ifndef ERF_SURFACE_LAYER_ROUGHNESS_H_
#define ERF_SURFACE_LAYER_ROUGHNESS_H_

#include <cmath>

#include "AMReX_GpuQualifiers.H"
#include "AMReX_REAL.H"

namespace surface_layer_roughness {

/**
 * Test whether a land surface model has supplied a usable momentum roughness
 * length for this cell.
 *
 * The LSM coupling MultiFabs are pre-filled with the `lsm_undefined` sentinel
 * (a very large positive number) and only overwritten where the provider
 * actually processed a column, so a sentinel is the normal state over water,
 * sea ice and any cell the provider skipped. On top of that, the MOST stress
 * kernels evaluate `log(zref/z0)` (ERF_MOSTStress.H), which is +inf at z0 == 0
 * and NaN for z0 < 0 -- so a non-positive roughness must never reach them.
 *
 * @param[in] znt        candidate roughness length from the LSM [m]
 * @param[in] undefined  the sentinel the LSM MultiFabs are filled with
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
lsm_roughness_is_valid(amrex::Real znt, amrex::Real undefined)
{
  return std::isfinite(znt) && (znt > amrex::Real(0.0)) &&
         (znt < amrex::Real(0.5) * undefined);
}

/**
 * Select the momentum roughness length for one surface cell.
 *
 * Over water the sea roughness models own z_0 and the LSM value is ignored.
 * Over land an LSM-provided value replaces the current one only when it is
 * usable; otherwise the existing z_0 (the constant, the value read from a
 * roughness file, or the value interpolated from level 0) is retained, so a
 * provider gap degrades to the pre-existing behavior rather than to a NaN.
 *
 * @param[in] is_land     true where the land mask marks this cell as land
 * @param[in] lsm_znt     candidate roughness length from the LSM [m]
 * @param[in] current_z0  the roughness length already stored for this cell [m]
 * @param[in] undefined   the sentinel the LSM MultiFabs are filled with
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
select_z0(
  bool is_land, amrex::Real lsm_znt, amrex::Real current_z0,
  amrex::Real undefined)
{
  if (!is_land) {
    return current_z0;
  }
  return lsm_roughness_is_valid(lsm_znt, undefined) ? lsm_znt : current_z0;
}

} // namespace surface_layer_roughness

#endif
