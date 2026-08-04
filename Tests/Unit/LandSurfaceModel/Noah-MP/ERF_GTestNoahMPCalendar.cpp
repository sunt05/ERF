/*
 * Unit tests for the Noah-MP calendar wiring (plan-slucm-coupling.md section 6,
 * Phase 0 prerequisite: "Wire JULIAN/YR/GMT").
 *
 * Before that wiring, NoahmpIOVarInitMod hard-set YR = 2000 / JULIAN = 1.0 and
 * nothing ever advanced them, so Noah-MP phenology under DYNAMIC_VEG_OPTION = 4
 * always interpolated January, the urban irrigation trigger never fired, and
 * SLUCM's shadow model would have had no clock.
 *
 * What is tested here is the exact arithmetic NOAHMP::Init and
 * NOAHMP::Advance_With_State use to stage YR / JULIAN / DECLIN:
 *   - erf_noahmp::epoch_from_calendar / calendar_from_epoch (defined in
 *     Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP_Init.cpp),
 *   - the orbital_decl() declination that Advance stages into NoahmpIO%DECLIN,
 *   - the solar hour angle OMG that NoahmpUrbanDriverMainMod derives from
 *     JULIAN + XLONG, checked against urban()'s tloc identity
 *     (module_sf_urban.F90:784 -- tloc == 12 at local solar noon).
 */

#include <cmath>
#include <string>

#include <gtest/gtest.h>

#include <AMReX.H>
#include <AMReX_BLassert.H>   // ERF_EpochTime.H uses AMREX_ASSERT but does not include this
#include <AMReX_REAL.H>

#include <ERF_Constants.H>
#include <ERF_EpochTime.H>
#include <ERF_OrbCosZenith.H>

// Declared the same way ERF_NOAHMP_Advance.cpp declares them (they have no
// header yet; see the note in ERF_NOAHMP_Init.cpp).
namespace erf_noahmp {
double epoch_from_calendar (int yr, double julian);
void   calendar_from_epoch (double epoch, int& yr, double& julian);
}

namespace
{

constexpr double sec_per_day = 86400.0;

// The WPS_Test configuration: namelist.erf START_* = 2023-08-17 06:00 UTC with
// NOAH_TIMESTEP = 3600.
constexpr int    wps_year   = 2023;
constexpr double wps_julian = 229.25;   // Aug 17 is day-of-year 229 (1-based); 06Z -> +0.25

// Exactly what NOAHMP::Advance_With_State computes each firing.
void stage_calendar (double epoch0, double elapsed_time,
                     int& yr, double& julian, double& declin)
{
    erf_noahmp::calendar_from_epoch(epoch0 + elapsed_time, yr, julian);

    double eccen = 0.0, obliq = 0.0, mvelp = 0.0;
    double obliqr = 0.0, lambm0 = 0.0, mvelpp = 0.0;
    int orbital_year = yr;
    orbital_params(orbital_year, eccen, obliq, mvelp, obliqr, lambm0, mvelpp);

    double eccf = 0.0, calday = julian;
    orbital_decl(calday, eccen, mvelpp, lambm0, obliqr, declin, eccf);
}

// NoahmpUrbanDriverMainMod's hour angle, and urban()'s tloc index built from it.
double hour_angle (double julian, double lon_deg)
{
    double omg = 2.0 * PI * (julian - std::floor(julian)) + lon_deg * PI / 180.0 - PI;
    // modulo(omg + pi, 2 pi) - pi
    double w = std::fmod(omg + PI, 2.0 * PI);
    if (w < 0.0) { w += 2.0 * PI; }
    return w - PI;
}

int urban_tloc (double omg)
{
    int t = static_cast<int>(omg / PI * 180.0 / 15.0 + 12.0 + 0.5) % 24;
    if (t < 0) { t += 24; }
    if (t == 0) { t = 24; }
    return t;
}

} // namespace

// The convention the whole chain depends on: JULIAN is 1-based day-of-year plus
// the UTC fraction of the day, so 1.0 is Jan 1 00Z. That is CAL_MON_DAY's
// convention (JULDAY = 1 -> Jan 1) and ERF_Radiation's `calday`.
TEST(NoahMPCalendar, JulianIsOneBasedDayOfYear)
{
    const double jan1 = static_cast<double>(getEpochTime("2023-01-01 00:00:00",
                                                         "%Y-%m-%d %H:%M:%S"));
    int yr = 0; double julian = 0.0;
    erf_noahmp::calendar_from_epoch(jan1, yr, julian);
    EXPECT_EQ(yr, 2023);
    EXPECT_NEAR(julian, 1.0, 1.0e-9);

    // 2023-08-17 06:00 UTC, the WPS_Test start date.
    const double aug17 = static_cast<double>(getEpochTime("2023-08-17 06:00:00",
                                                          "%Y-%m-%d %H:%M:%S"));
    erf_noahmp::calendar_from_epoch(aug17, yr, julian);
    EXPECT_EQ(yr, wps_year);
    EXPECT_NEAR(julian, wps_julian, 1.0e-9);
}

// epoch_from_calendar is the inverse used by Advance to rebuild the origin epoch
// from the START_YR / START_JULIAN pinned into every NoahmpIO block at Init.
TEST(NoahMPCalendar, EpochRoundTrip)
{
    for (int yr : {1999, 2000, 2020, 2023, 2024, 2100}) {
        for (double j : {1.0, 59.5, 60.0, 180.125, 365.75}) {
            const double epoch = erf_noahmp::epoch_from_calendar(yr, j);
            int    yr_rt = 0;
            double j_rt  = 0.0;
            erf_noahmp::calendar_from_epoch(epoch, yr_rt, j_rt);
            EXPECT_EQ(yr_rt, yr) << " julian " << j;
            EXPECT_NEAR(j_rt, j, 1.0e-6) << " year " << yr;
        }
    }
}

// The regression this whole change exists for: JULIAN must advance, one
// NOAH_TIMESTEP at a time, instead of sitting at 1.0 forever.
TEST(NoahMPCalendar, JulianAdvancesEveryTimestep)
{
    const double epoch0 = erf_noahmp::epoch_from_calendar(wps_year, wps_julian);
    const double dtbl   = 3600.0;   // namelist.erf NOAH_TIMESTEP for WPS_Test

    double prev_julian = 0.0;
    for (int step = 0; step <= 48; ++step) {
        int    yr     = 0;
        double julian = 0.0, declin = 0.0;
        stage_calendar(epoch0, static_cast<double>(step) * dtbl, yr, julian, declin);

        EXPECT_EQ(yr, wps_year);
        EXPECT_NEAR(julian, wps_julian + step * dtbl / sec_per_day, 1.0e-9);
        if (step > 0) {
            EXPECT_GT(julian, prev_julian);
            EXPECT_NEAR(julian - prev_julian, dtbl / sec_per_day, 1.0e-9);
        }
        prev_julian = julian;

        // Mid-August declination: ~ +13.5 deg, and falling through the period.
        EXPECT_GT(declin, 0.20);
        EXPECT_LT(declin, 0.26);
    }

    // Two days on, the phenology month index has moved off the origin: this is
    // what DYNAMIC_VEG_OPTION = 4 keys on (PhenologyMainMod interpolates the
    // monthly LAI table from DayJulianInYear).
    int    yr_end = 0;
    double julian_end = 0.0, declin_end = 0.0;
    stage_calendar(epoch0, 2.0 * sec_per_day, yr_end, julian_end, declin_end);
    EXPECT_NEAR(julian_end, wps_julian + 2.0, 1.0e-9);
    EXPECT_GT(12.0 * julian_end / 365.0, 7.0);   // month index ~7.5, i.e. August
    EXPECT_LT(12.0 * julian_end / 365.0, 8.0);
}

// A run that crosses New Year must roll the year over and restart day-of-year at
// 1.0; a leap year must have 366 days.
TEST(NoahMPCalendar, YearRolloverAndLeapYear)
{
    // 2023-12-31 23:00 UTC + 2 h
    const double epoch0 = erf_noahmp::epoch_from_calendar(2023, 365.0 + 23.0 / 24.0);
    int    yr = 0;
    double julian = 0.0, declin = 0.0;

    stage_calendar(epoch0, 0.0, yr, julian, declin);
    EXPECT_EQ(yr, 2023);
    EXPECT_NEAR(julian, 365.0 + 23.0 / 24.0, 1.0e-9);

    stage_calendar(epoch0, 7200.0, yr, julian, declin);
    EXPECT_EQ(yr, 2024);
    EXPECT_NEAR(julian, 1.0 + 1.0 / 24.0, 1.0e-9);

    // 2024 is a leap year: Feb 29 exists and Dec 31 is day 366.
    const double leap_epoch = erf_noahmp::epoch_from_calendar(2024, 1.0);
    stage_calendar(leap_epoch, 365.0 * sec_per_day, yr, julian, declin);
    EXPECT_EQ(yr, 2024);
    EXPECT_NEAR(julian, 366.0, 1.0e-9);
}

// The solar hour angle NoahmpUrbanDriverMainMod derives for SLUCM. urban() turns
// it back into an hour index with
//   tloc = mod(int(OMG/PI*180./15.+12.+0.5),24)   (module_sf_urban.F90:784)
// which must read 12 at local solar noon for the anthropogenic-heat and
// irrigation diurnal profiles to be indexed correctly.
TEST(NoahMPCalendar, HourAngleIsTwelveAtLocalSolarNoon)
{
    // Greenwich: local solar noon is 12Z.
    EXPECT_NEAR(hour_angle(200.5, 0.0), 0.0, 1.0e-12);
    EXPECT_EQ(urban_tloc(hour_angle(200.5, 0.0)), 12);

    // Every whole-degree longitude: local solar noon is at UTC = 12 - lon/15.
    for (int lon_deg = -180; lon_deg <= 180; lon_deg += 5) {
        const double utc_noon = 12.0 - lon_deg / 15.0;      // may fall outside [0,24)
        const double julian   = 200.0 + utc_noon / 24.0;
        const double omg      = hour_angle(julian, static_cast<double>(lon_deg));
        EXPECT_NEAR(omg, 0.0, 1.0e-9) << " lon " << lon_deg;
        EXPECT_EQ(urban_tloc(omg), 12) << " lon " << lon_deg;
    }

    // ... and the index tracks local solar time hour by hour.
    for (int h = 0; h < 24; ++h) {
        const double julian = 200.0 + h / 24.0;             // UTC hour h
        const int    tloc   = urban_tloc(hour_angle(julian, 0.0));
        EXPECT_EQ(tloc, (h == 0) ? 24 : h) << " utc hour " << h;
    }

    // OMG is normalised to (-pi, pi] whatever the longitude/day fraction.
    for (int lon_deg = -180; lon_deg <= 180; lon_deg += 15) {
        for (int h = 0; h < 24; ++h) {
            const double omg = hour_angle(200.0 + h / 24.0, static_cast<double>(lon_deg));
            EXPECT_GT(omg, -PI - 1.0e-12);
            EXPECT_LE(omg,  PI + 1.0e-12);
        }
    }
}

// Declination must follow the season, not sit at whatever Jan 1 gives. This is
// the property the frozen calendar destroyed for SLUCM's shadow geometry.
TEST(NoahMPCalendar, DeclinationFollowsTheSeason)
{
    auto declin_on = [](int yr, double julian) {
        int    y = 0;
        double j = 0.0, d = 0.0;
        stage_calendar(erf_noahmp::epoch_from_calendar(yr, julian), 0.0, y, j, d);
        return d;
    };

    const double d_jan = declin_on(2023,   1.5);   // Jan 1  12Z
    const double d_mar = declin_on(2023,  80.5);   // ~equinox
    const double d_jun = declin_on(2023, 172.5);   // ~solstice
    const double d_dec = declin_on(2023, 355.5);   // ~solstice

    EXPECT_LT(d_jan, -0.36);  EXPECT_GT(d_jan, -0.41);   // ~ -23 deg
    EXPECT_NEAR(d_mar, 0.0, 0.02);                        // ~ 0 deg
    EXPECT_GT(d_jun,  0.39);  EXPECT_LT(d_jun,  0.42);   // ~ +23.4 deg
    EXPECT_LT(d_dec, -0.39);  EXPECT_GT(d_dec, -0.42);   // ~ -23.4 deg

    // The frozen-calendar bug in one assertion: staging Jan 1 where the model is
    // really in June is a ~47 degree declination error.
    EXPECT_GT(std::abs(d_jun - d_jan) * 180.0 / PI, 45.0);
}
