"""
psefinder.py
------------
Searches for potential solar eclipses within a specified time interval
using Skyfield ephemerides and custom project utilities.
"""

import math

from skyfield.api import load, GREGORIAN_START

import psenarrow       # Precise eclipse timing
import pconstants      # Physical constants like MOON_RADIUS_KM, SUN_RADIUS_KM
import pdefilepath     # Ephemeris file paths
import pedatetime      # Custom datetime class


def sefinder(
    start_time: pedatetime.datetime,
    end_time: pedatetime.datetime,
    step: pedatetime.timedelta,
    tt_enable: bool = False,
    printsep: bool = False,
) -> None:
    """
    Search for potential solar eclipses within a given time range.

    Parameters
    ----------
    start_time : pedatetime.datetime
        Start of the search interval.
    end_time : pedatetime.datetime
        End of the search interval.
    step : pedatetime.timedelta
        Step size for scanning the interval.
    tt_enable : bool, default False
        Whether to apply TT/ΔT correction to the refined eclipse time.
    printsep : bool, default False
        Whether to print the minimum angular separation along with the date.

    Notes
    -----
    The function scans the interval in steps of `step` and detects potential
    eclipses based on angular separation and geometric thresholds. Detected
    eclipses are refined using `psenarrow.senarrow`.
    """
    # Load ephemerides and timescale
    eph = load(pdefilepath.EPHEM_PATH)
    ts = load.timescale()
    ts.julian_calendar_cutoff = GREGORIAN_START

    earth, sun, moon = eph["earth"], eph["sun"], eph["moon"]

    current_time = start_time.copy()

    while current_time <= end_time:
        # Convert current time to Skyfield Time object (UT1)
        sf_time = ts.ut1(
            current_time.year,
            current_time.month,
            current_time.day,
            current_time.hour,
            current_time.minute,
            current_time.second,
        )

        if tt_enable:
            sf_time = ts.tt(
                current_time.year,
                current_time.month,
                current_time.day,
                current_time.hour,
                current_time.minute,
                current_time.second,
            )

        # Apparent positions of Sun and Moon
        sun_pos = earth.at(sf_time).observe(sun).apparent()
        moon_pos = earth.at(sf_time).observe(moon).apparent()

        # Angular separation in radians
        sep_angle: float = moon_pos.separation_from(sun_pos).radians

        # Distances in kilometers
        sun_dist_km: float = sun_pos.distance().km
        moon_dist_km: float = moon_pos.distance().km

        # Eclipse geometry threshold in radians
        threshold: float = (
            math.asin((pconstants.MOON_RADIUS_KM + pconstants.EARTH_RADIUS_KM) / moon_dist_km)
            + math.asin((pconstants.SUN_RADIUS_KM - pconstants.EARTH_RADIUS_KM) / sun_dist_km)
        )

        # Potential eclipse detected
        if sep_angle <= threshold:
            # Estimate start slightly before current scan
            start_est = current_time.copy()
            start_est.sub_minute()

            # Estimate end a few hours after current scan
            end_est = current_time.copy()
            end_est.add_hours(3)

            # Refine eclipse time and minimum separation
            eclipse_date, min_sep = psenarrow.senarrow(start_est, end_est, tt_enable)

            # Add time system suffix for clarity
            if eclipse_date is not None:
                suffix: str = " TT" if tt_enable else " UT1"
                eclipse_date = eclipse_date + suffix

                # Print results
                if printsep:
                    print(f"{eclipse_date}, {min_sep} rad")
                else:
                    print(f"{eclipse_date}")

            # Skip ~27 days to avoid detecting the same eclipse again
            current_time.add_days(27)
        else:
            # Advance by the specified step
            current_time = current_time + step
