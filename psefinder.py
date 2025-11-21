"""
psefinder.py
-----------------

Searches for potential solar eclipses within a specified time interval using
Skyfield ephemerides and custom project utilities.

"""

from skyfield.api import load
import math

import psenarrow
import pconstants
import pdefilepath
import pedatetime


def sefinder(
    start_time: pedatetime.datetime, end_time: pedatetime.datetime, printsep: bool
) -> None:
    """
    Search for possible solar eclipses within a given time range.

    Parameters
    ----------
    start_time : pedatetime.datetime
        The start of the search interval.
    end_time : pedatetime.datetime
        The end of the search interval.
    printsep : bool
        Whether to print the minimum angular separation along with the date.

    Notes
    -----
    This function scans the time range in 2-hour steps to detect potential
    eclipses based on the angular separation between the Sun and Moon and
    a geometric threshold.
    """

    # Load planetary ephemeris from file (e.g., DE441)
    planets = load(pdefilepath.EPHEM_PATH)
    # Load timescale object for date/time conversions
    ts = load.timescale()

    # Extract Earth, Sun, and Moon objects
    earth, sun, moon = planets["earth"], planets["sun"], planets["moon"]

    # Initialize current scanning time
    curtime = start_time.copy()

    while curtime <= end_time:
        # Convert current time to Skyfield Time object
        t = ts.utc(
            curtime.year,
            curtime.month,
            curtime.day,
            curtime.hour,
            curtime.minute,
            curtime.second,
        )

        # Compute apparent positions of Sun and Moon from Earth
        obs_sun = earth.at(t).observe(sun)
        obs_moon = earth.at(t).observe(moon)

        # Compute angular separation between Sun and Moon (in radians)
        sep_angle = obs_moon.apparent().separation_from(obs_sun.apparent()).radians

        # Get distances to Sun and Moon in kilometers
        _, _, sdist = obs_sun.radec()
        _, _, mdist = obs_moon.radec()
        sdist_km, mdist_km = sdist.km, mdist.km

        # Compute the eclipse geometry threshold (radians)
        # This accounts for the apparent radii of Sun, Moon, and Earth
        eclipse_threshold_angle = math.asin(
            (pconstants.MOON_RADIUS_KM + pconstants.EARTH_RADIUS_KM) / mdist_km
        ) + math.asin(
            (pconstants.SUN_RADIUS_KM - pconstants.EARTH_RADIUS_KM) / sdist_km
        )

        # Check if angular separation is less than threshold (possible eclipse)
        # If so, narrow down to the exact start and end times
        if eclipse_threshold_angle >= sep_angle:
            # Estimate start time of eclipse slightly earlier
            tstart = curtime.copy()
            tstart.sub_minute()

            # Estimate end time of eclipse a few hours later
            tend = curtime.copy()
            tend.add_hours(3)

            # Compute the precise time and minimum separation using senarrow
            date, min_sep = psenarrow.senarrow(tstart, tend)

            if date:
                if not printsep:
                    # Only print the date of maximum eclipse
                    print(date)
                else:
                    # Print date and minimum angular separation
                    print(date, f"{min_sep} rad")

            # Skip ahead roughly one synodic month (~27 days)
            # to avoid redundant scanning of the same eclipse
            curtime.add_days(27)

        # Step forward 2 hours (smaller steps = more accurate, slower)
        curtime.add_hours(2)
