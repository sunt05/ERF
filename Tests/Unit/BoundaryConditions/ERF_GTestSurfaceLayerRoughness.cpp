#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "AMReX_REAL.H"

#include "ERF_Constants.H"
#include "ERF_SurfaceLayerRoughness.H"

using amrex::Real;

namespace {

constexpr Real undefined = lsm_undefined;

// A plausible pre-existing MOST roughness (erf.most.z0 default) and a
// plausible Noah-MP land value, deliberately different so a silent no-op is
// distinguishable from a successful update.
constexpr Real current_z0 = Real(0.1);
constexpr Real lsm_z0 = Real(0.47);

constexpr bool land = true;
constexpr bool water = false;

} // namespace

// Motivation: the whole point of roughness_type_land = lsm is that the land
// surface model's per-step roughness reaches MOST; over land a usable value
// must replace whatever z_0 currently holds.
TEST(SurfaceLayerRoughness, ValidLandValueReplacesCurrentRoughness)
{
  EXPECT_TRUE(
    surface_layer_roughness::lsm_roughness_is_valid(lsm_z0, undefined));
  EXPECT_EQ(
    surface_layer_roughness::select_z0(land, lsm_z0, current_z0, undefined),
    lsm_z0);
}

// Motivation: the LSM coupling MultiFabs are pre-filled with lsm_undefined and
// only overwritten where the provider ran, so the sentinel is the normal state
// before the first LSM call and over every unprocessed column. Accepting it
// would drive log(zref/z0) to a large negative value and invert the stress.
TEST(SurfaceLayerRoughness, UndefinedSentinelIsRejectedAndKeepsCurrentRoughness)
{
  EXPECT_FALSE(
    surface_layer_roughness::lsm_roughness_is_valid(undefined, undefined));
  EXPECT_EQ(
    surface_layer_roughness::select_z0(land, undefined, current_z0, undefined),
    current_z0);
}

// Motivation: result_is_valid in the Noah-MP policy accepts 0.0, but
// ERF_MOSTStress.H evaluates log(zref/z0), which is +inf at z0 == 0 and NaN
// for z0 < 0. The roughness selection must reject both itself.
TEST(SurfaceLayerRoughness, NonPositiveRoughnessIsRejected)
{
  for (const Real bad : {Real(0.0), Real(-0.0), Real(-0.1), Real(-9999.0)}) {
    EXPECT_FALSE(surface_layer_roughness::lsm_roughness_is_valid(bad, undefined))
      << "value " << bad;
    EXPECT_EQ(
      surface_layer_roughness::select_z0(land, bad, current_z0, undefined),
      current_z0)
      << "value " << bad;
  }
}

// Motivation: a NaN or infinite provider value would propagate silently through
// log(zref/z0) into u* and every downstream flux.
TEST(SurfaceLayerRoughness, NonFiniteRoughnessIsRejected)
{
  for (const Real bad :
       {std::numeric_limits<Real>::infinity(),
        -std::numeric_limits<Real>::infinity(),
        std::numeric_limits<Real>::quiet_NaN()}) {
    EXPECT_FALSE(surface_layer_roughness::lsm_roughness_is_valid(bad, undefined));
    EXPECT_EQ(
      surface_layer_roughness::select_z0(land, bad, current_z0, undefined),
      current_z0);
  }
}

// Motivation: over water the Charnock/Donelan/wave-coupled models own z_0 and
// update it inside the flux iteration. An LSM value must never overwrite it,
// even when the value itself is perfectly usable.
TEST(SurfaceLayerRoughness, WaterCellsKeepTheirSeaRoughness)
{
  EXPECT_EQ(
    surface_layer_roughness::select_z0(water, lsm_z0, current_z0, undefined),
    current_z0);
  EXPECT_EQ(
    surface_layer_roughness::select_z0(water, undefined, current_z0, undefined),
    current_z0);
}

// Motivation: the sentinel test is a magnitude cutoff, not an equality test, so
// that partially-overwritten or scaled sentinels are also caught. Physically
// plausible roughness lengths must stay comfortably inside the accepted range.
TEST(SurfaceLayerRoughness, CutoffAcceptsPhysicalRoughnessAndRejectsHugeValues)
{
  // Smooth water through dense urban canopy.
  for (const Real good :
       {Real(1.0e-5), Real(0.01), Real(0.1), Real(0.5), Real(2.0)}) {
    EXPECT_TRUE(surface_layer_roughness::lsm_roughness_is_valid(good, undefined))
      << "value " << good;
  }
  EXPECT_FALSE(surface_layer_roughness::lsm_roughness_is_valid(
    Real(0.6) * undefined, undefined));
  EXPECT_TRUE(surface_layer_roughness::lsm_roughness_is_valid(
    Real(0.4) * undefined, undefined));
}

// Motivation: z_0 is updated in place every step, so the selection is applied
// repeatedly to its own output. A provider gap after a good step must hold the
// last accepted LSM value rather than snapping back to the namelist constant.
TEST(SurfaceLayerRoughness, SelectionIsStableAcrossRepeatedApplication)
{
  const Real after_good =
    surface_layer_roughness::select_z0(land, lsm_z0, current_z0, undefined);
  EXPECT_EQ(after_good, lsm_z0);

  const Real after_gap =
    surface_layer_roughness::select_z0(land, undefined, after_good, undefined);
  EXPECT_EQ(after_gap, lsm_z0);

  const Real after_good_again =
    surface_layer_roughness::select_z0(land, lsm_z0, after_gap, undefined);
  EXPECT_EQ(after_good_again, lsm_z0);
}
