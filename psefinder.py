"""
psefinder.py
-----------------

Searches for potential solar eclipses within a specified time interval using
Skyfield ephemerides and custom project utilities.
"""

from skyfield.api import load
import math

import psenarrow  # Custom module for precise eclipse timing
import pconstants  # Custom module with constants like MOON_RADIUS_KM, SUN_RADIUS_KM
import pdefilepath  # Custom module with ephemeris file paths
import pedatetime  # Custom datetime class used in this project


def sefinder(
    start_time: pedatetime.datetime,
    end_time: pedatetime.datetime,
    printsep: bool,
    step: pedatetime.timedelta,
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
    shrs : int
        Step size in hours for scanning the time interval.

    Notes
    -----
    This function scans the time range in `shrs`-hour steps to detect potential
    eclipses based on the angular separation between the Sun and Moon and
    a geometric threshold.
    """

    # Load planetary ephemeris (e.g., DE441)
    eph = load(pdefilepath.EPHEM_PATH)

    # Load Skyfield timescale for converting datetimes
    timescale = load.timescale()

    # Extract Earth, Sun, and Moon objects for position calculations
    earth, sun, moon = eph["earth"], eph["sun"], eph["moon"]

    # Initialize current scanning time
    current_time = start_time.copy()

    # Loop through time interval
    while current_time <= end_time:
        # Convert current time to Skyfield Time object
        skyfield_time = timescale.utc(
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

        # Compute angular separation (in radians) between Sun and Moon
        angular_separation = (
            moon_position.apparent().separation_from(sun_position.apparent()).radians
        )

        # Retrieve distances to Sun and Moon (in km)
        _, _, sun_distance = sun_position.radec()
        _, _, moon_distance = moon_position.radec()
        sun_distance_km, moon_distance_km = sun_distance.km, moon_distance.km

        # Compute eclipse geometry threshold (radians)
        # Based on apparent radii of Moon and Sun and Earth's radius
        eclipse_threshold = math.asin(
            (pconstants.MOON_RADIUS_KM + pconstants.EARTH_RADIUS_KM) / moon_distance_km
        ) + math.asin(
            (pconstants.SUN_RADIUS_KM - pconstants.EARTH_RADIUS_KM) / sun_distance_km
        )

        # Check if current separation indicates a possible eclipse
        if eclipse_threshold >= angular_separation:
            # Estimate start time slightly before current scan time
            start_estimate = current_time.copy()
            start_estimate.sub_minute()

            # Estimate end time a few hours after current scan time
            end_estimate = current_time.copy()
            end_estimate.add_hours(3)

            # Refine to precise eclipse time and minimum angular separation
            eclipse_date, min_separation = psenarrow.senarrow(
                start_estimate, end_estimate
            )

            # Print results if a valid eclipse is found
            if eclipse_date:
                if not printsep:
                    # Only print date of maximum eclipse
                    print(eclipse_date)
                else:
                    # Print date and minimum angular separation
                    print(eclipse_date, f"{min_separation} rad")

            # Skip roughly one synodic month (~27 days) to avoid detecting same eclipse again
            current_time.add_days(27)

        # Step forward by the specified timedelta
        current_time = current_time + step
