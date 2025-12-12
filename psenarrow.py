"""
psenarrow.py
-----------------

Provides a high-resolution search routine for locating the exact moment of
closest angular approach between the Sun and Moon.
"""

from skyfield.api import load, GREGORIAN_START
import math

import pconstants  # Custom module containing constants like MOON_RADIUS_KM, EARTH_RADIUS_KM, SUN_RADIUS_KM
import pdefilepath  # Custom module containing file paths for ephemerides
import pedatetime  # Custom datetime class used in this project


def senarrow(
    starttime: pedatetime.datetime, endtime: pedatetime.datetime, tt: bool
) -> tuple[str, float]:
    """
    Finds the time and angular distance of the closest Sun-Moon approach
    within the given time interval.

    Args:
        starttime (pedatetime.datetime): Start of search range (UTC)
        endtime (pedatetime.datetime): End of search range (UTC)

    Returns:
        tuple:
            date (str): Time of minimum separation in ISO format
            min_sep_angle (float): Minimum angular separation (radians)
            Returns (None, None) if no eclipse is detected.
    """

    # Lists to store angular separations and corresponding timestamps
    angular_separations = []
    timestamps = []

    # Load planetary ephemerides and timescale
    eph = load(pdefilepath.EPHEM_PATH)
    ts = load.timescale()
    ts.julian_calendar_cutoff = GREGORIAN_START

    # Extract Earth, Sun, and Moon objects
    earth, sun, moon = eph["earth"], eph["sun"], eph["moon"]

    # Initialize current scanning time
    current_time = starttime.copy()

    # Iterate through each second in the interval
    while current_time <= endtime:
        # Convert current time to Skyfield Time object
        skyfield_time = ts.utc(
            current_time.year,
            current_time.month,
            current_time.day,
            current_time.hour,
            current_time.minute,
            current_time.second,
        )

        # Compute apparent positions of Sun and Moon from Earth
        sun_position = earth.at(skyfield_time).observe(sun)
        moon_position = earth.at(skyfield_time).observe(moon)

        # Compute angular separation (radians) between Moon and Sun
        angular_separation = (
            moon_position.apparent().separation_from(sun_position.apparent()).radians
        )

        # Get distances of Sun and Moon from Earth in kilometers
        _, _, sun_distance = sun_position.radec()
        _, _, moon_distance = moon_position.radec()
        sun_distance_km, moon_distance_km = sun_distance.km, moon_distance.km

        # Compute eclipse geometry threshold (radians)
        eclipse_threshold = math.asin(
            (pconstants.MOON_RADIUS_KM + pconstants.EARTH_RADIUS_KM) / moon_distance_km
        ) + math.asin(
            (pconstants.SUN_RADIUS_KM - pconstants.EARTH_RADIUS_KM) / sun_distance_km
        )

        # Only consider times where separation is within eclipse threshold
        if eclipse_threshold >= angular_separation:
            angular_separations.append(angular_separation)  # Store angular separation
            timestamps.append(current_time.copy())  # Store timestamp

        # Move to the next second
        current_time.add_second()

    # Edge case: No times were within threshold
    if not angular_separations:
        return None, None

    # Find the minimum separation and corresponding timestamp
    min_sep_angle = min(angular_separations)
    min_index = angular_separations.index(min_sep_angle)
    min_time = timestamps[min_index]

    # Add tt if the user wants
    t = ts.utc(
        min_time.year,
        min_time.month,
        min_time.day,
        min_time.hour,
        min_time.minute,
        min_time.second,
    )

    return_time = min_time.copy()
    if tt:
        return_time = min_time.copy() + pedatetime.timedelta(
            0, 0, 0, round(float(t.delta_t))
        )
    else:
        return_time = min_time.copy()

    # Return ISO format time and minimum separation
    return return_time.isoformat(), min_sep_angle
