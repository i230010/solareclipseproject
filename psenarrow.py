"""
psenarrow.py
------------
High-resolution search for the exact moment of closest angular approach
between the Sun and Moon within a given time interval.
"""

import math
from typing import List, Optional, Tuple

from skyfield.api import load, GREGORIAN_START

import pconstants  # Physical constants (MOON_RADIUS_KM, SUN_RADIUS_KM, EARTH_RADIUS_KM)
import pdefilepath  # Ephemeris file path
import pedatetime  # Custom datetime class


def senarrow(
    starttime: pedatetime.datetime,
    endtime: pedatetime.datetime,
    tt_enable: bool = False,
) -> Tuple[Optional[str], Optional[float]]:
    """
    Find the time and angular distance of the closest Sun-Moon approach
    (potential solar eclipse) within a given interval.

    Parameters
    ----------
    starttime : pedatetime.datetime
        Start of the search interval (UTC)
    endtime : pedatetime.datetime
        End of the search interval (UTC)
    tt_enable : bool, default False
        Whether to apply TT/delta T correction to the resulting time.

    Returns
    -------
    Tuple[Optional[str], Optional[float]]
        ISO-format time of minimum separation and minimum angular separation
        in radians. Returns (None, None) if no eclipse is detected.
    """
    if starttime > endtime:
        raise ValueError("starttime must be earlier than or equal to endtime")

    # Load ephemerides and timescale
    eph = load(pdefilepath.EPHEM_PATH)
    ts = load.timescale()
    ts.julian_calendar_cutoff = GREGORIAN_START

    earth, sun, moon = eph["earth"], eph["sun"], eph["moon"]

    separations: List[float] = []
    timestamps: List[pedatetime.datetime] = []

    current_time = starttime.copy()

    # Scan each second in the interval for closest approach
    while current_time <= endtime:
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

        # Apparent positions
        sun_pos = earth.at(sf_time).observe(sun).apparent()
        moon_pos = earth.at(sf_time).observe(moon).apparent()

        # Angular separation in radians
        sep_angle: float = moon_pos.separation_from(sun_pos).radians

        # Distances to Sun and Moon in kilometers
        sun_dist_km: float = sun_pos.distance().km
        moon_dist_km: float = moon_pos.distance().km

        # Eclipse threshold in radians based on apparent sizes
        threshold: float = math.asin(
            (pconstants.MOON_RADIUS_KM + pconstants.EARTH_RADIUS_KM) / moon_dist_km
        ) + math.asin(
            (pconstants.SUN_RADIUS_KM - pconstants.EARTH_RADIUS_KM) / sun_dist_km
        )

        if sep_angle <= threshold:
            separations.append(sep_angle)
            timestamps.append(current_time.copy())

        # Move forward by 1 second
        current_time.add_second()

    # If no eclipse detected
    if not separations:
        return None, None

    # Identify minimum angular separation and corresponding time
    min_sep = min(separations)
    min_index = separations.index(min_sep)
    min_time = timestamps[min_index]

    return min_time.isoformat(), min_sep
