#!/usr/bin/env python3
"""Generate the 48 h synthetic forcing for the SLUCM single-column harness.

Every field urban() reads as INTENT(IN) is written explicitly, so the ERF and
upstream builds cannot diverge through the driver. Values are emitted with 17
significant digits, which round-trips an IEEE double exactly -- both builds read
identical bits.

The column is a mid-latitude summer day: a smooth diurnal temperature wave, a
clear-sky shortwave cycle from the solar geometry, constant humidity, wind and
longwave, and no rain. It is not a real site; item 5 of the plan's section 8
(the Grimmond et al. urban energy-balance comparison) is what brings observed
forcing, and this file's column layout is the one that reader will produce.

    ./make_forcing.py > forcing_48h.csv
"""

import math
import sys

DELT = 60.0          # s
HOURS = 48
NSTEP = int(HOURS * 3600 / DELT)

LAT_DEG = 40.0       # must match XLAT_ in slucm_column.F90
DECLIN = 0.4         # rad, ~summer
LON_DEG = 0.0

TA_MEAN = 291.0      # K
TA_AMP = 6.0         # K
QA = 0.008           # kg/kg
UA = 3.0             # m/s
LLG = 320.0          # W/m2, downward longwave
RHOO = 1.2           # kg/m3
S0 = 900.0           # W/m2, clear-sky peak normal irradiance


def fmt(x):
    return repr(float(x)) if False else "%.17e" % x


def main():
    lat = math.radians(LAT_DEG)
    out = sys.stdout
    out.write("DELT,TA,QA,UA,U1,V1,SSG,LLG,RAIN,RHOO,COSZ,OMG,DECLIN\n")

    for n in range(NSTEP):
        t = (n + 1) * DELT                 # s since start
        hour = (t / 3600.0) % 24.0         # local solar hour

        # hour angle: 0 at local solar noon, negative before
        omg = math.radians(15.0 * (hour - 12.0)) + math.radians(LON_DEG)
        omg = (omg + math.pi) % (2.0 * math.pi) - math.pi

        cosz = math.sin(lat) * math.sin(DECLIN) + \
            math.cos(lat) * math.cos(DECLIN) * math.cos(omg)
        cosz = max(cosz, 0.0)

        ssg = S0 * cosz                    # downward shortwave at the surface
        ta = TA_MEAN + TA_AMP * math.sin(2.0 * math.pi * (hour - 9.0) / 24.0)

        row = [DELT, ta, QA, UA, UA, 0.0, ssg, LLG, 0.0, RHOO, cosz, omg, DECLIN]
        out.write(",".join(fmt(v) for v in row) + "\n")


if __name__ == "__main__":
    main()
