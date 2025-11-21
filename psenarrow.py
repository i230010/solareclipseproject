"""
psenarrow.py
-----------------

Provides a high-resolution search routine for locating the exact moment of
closest angular approach between the Sun and Moon.

"""

from skyfield.api import load
import math

import pconstants
import pdefilepath
import pedatetime


def senarrow(starttime: pedatetime.datetime, endtime: pedatetime.datetime):
    """
    Finds the time and angular distance of the closest Sun-Moon approach.

    Args:
            starttime (datetime): Start of search range (UTC)
            endtime (datetime): End of search range (UTC)

    Returns:
            tuple:
                    date (datetime): Time of minimum separation
                    min_sep_angle (float): Minimum angular separation (radians)
    """
    separations = []
    timestamps = []

    # Load Skyfield ephemerides and timescale once
    planets = load(pdefilepath.EPHEM_PATH)
    ts = load.timescale()

    earth, sun, moon = planets["earth"], planets["sun"], planets["moon"]

    curtime = starttime.copy()

    # Step through each second to find region of closest approach
    while curtime <= endtime:
        t = ts.utc(
            curtime.year,
            curtime.month,
            curtime.day,
            curtime.hour,
            curtime.minute,
            curtime.second,
        )

        # Apparent positions
        obs_sun = earth.at(t).observe(sun)
        obs_moon = earth.at(t).observe(moon)

        # Angular separation (radians)
        sep_angle = obs_moon.apparent().separation_from(obs_sun.apparent()).radians

        # Distances (km)
        _, _, sdist = obs_sun.radec()
        _, _, mdist = obs_moon.radec()
        sdist_km, mdist_km = sdist.km, mdist.km

        # Eclipse geometry threshold (radians)
        eclipse_threshold_angle = math.asin(
            (pconstants.MOON_RADIUS_KM + pconstants.EARTH_RADIUS_KM) / mdist_km
        ) + math.asin(
            (pconstants.SUN_RADIUS_KM - pconstants.EARTH_RADIUS_KM) / sdist_km
        )

        if eclipse_threshold_angle >= sep_angle:
            separations.append(sep_angle)
            timestamps.append(curtime.copy())

        # Increment time by 1 second
        curtime.add_second()

    min_sep_angle = min(separations)
    min_index = separations.index(min_sep_angle)
    date = timestamps[min_index]

    return date.isoformat(), min_sep_angle
