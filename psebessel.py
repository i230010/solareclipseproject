"""
psebessel.py
-------------

Computes Besselian elements for a given datetime and fits
polynomials to these elements for smooth interpolation.
"""

import math
import numpy as np
from skyfield.api import load, GREGORIAN_START
from skyfield.units import Angle
import vector

import pconstants
import pdefilepath
import pedatetime


# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------


def km_to_earth_radii(km: float) -> float:
    """
    Convert a distance in kilometers to Earth radii.
    """
    return km / pconstants.EARTH_RADIUS_KM


# ---------------------------------------------------------------------------
# Compute Besselian Elements
# ---------------------------------------------------------------------------


def besselian_find(
    dt: pedatetime.datetime,
) -> tuple[float, float, float, float, float, float, float, float]:
    """
    Compute Besselian elements for a given datetime.

    Returns 8 floats:
        moon_x, moon_y, shadow_decl_deg, northern_limit, southern_limit,
        sun_hour_angle_deg, tangent_north, tangent_south
    """
    planets = load(pdefilepath.EPHEM_PATH)
    earth, sun, moon = planets["earth"], planets["sun"], planets["moon"]

    ts = load.timescale()
    ts.julian_calendar_cutoff = GREGORIAN_START
    t_sf = ts.tt(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)

    obs_sun = earth.at(t_sf).observe(sun)
    obs_moon = earth.at(t_sf).observe(moon)

    sun_ra, sun_dec, sun_dist = obs_sun.radec()
    moon_ra, moon_dec, moon_dist = obs_moon.radec()

    sun_radius_r = km_to_earth_radii(sun_dist.km)
    moon_radius_r = km_to_earth_radii(moon_dist.km)

    sun_ra_rad, sun_dec_rad = sun_ra.radians, sun_dec.radians
    moon_ra_rad, moon_dec_rad = moon_ra.radians, moon_dec.radians

    # Rectangular coordinates
    sun_vec = vector.obj(
        x=sun_radius_r * math.cos(sun_dec_rad) * math.cos(sun_ra_rad),
        y=sun_radius_r * math.cos(sun_dec_rad) * math.sin(sun_ra_rad),
        z=sun_radius_r * math.sin(sun_dec_rad),
    )
    moon_vec = vector.obj(
        x=moon_radius_r * math.cos(moon_dec_rad) * math.cos(moon_ra_rad),
        y=moon_radius_r * math.cos(moon_dec_rad) * math.sin(moon_ra_rad),
        z=moon_radius_r * math.sin(moon_dec_rad),
    )

    shadow_vec = sun_vec - moon_vec
    shadow_dist = abs(shadow_vec)
    shadow_axis_angle = math.atan2(shadow_vec.y, shadow_vec.x)
    shadow_decl = math.asin(shadow_vec.z / shadow_dist)
    sun_hour_angle = (Angle(degrees=t_sf.gmst * 15).radians - shadow_axis_angle) % (
        2.0 * math.pi
    )

    moon_x = moon_radius_r * (
        math.cos(moon_dec_rad) * math.sin(moon_ra_rad - shadow_axis_angle)
    )
    moon_y = moon_radius_r * (
        math.sin(moon_dec_rad) * math.cos(shadow_decl)
        - math.cos(moon_dec_rad)
        * math.sin(shadow_decl)
        * math.cos(moon_ra_rad - shadow_axis_angle)
    )
    moon_z = moon_radius_r * (
        math.sin(moon_dec_rad) * math.sin(shadow_decl)
        + math.cos(moon_dec_rad)
        * math.cos(shadow_decl)
        * math.cos(moon_ra_rad - shadow_axis_angle)
    )

    sun_radius = km_to_earth_radii(pconstants.SUN_RADIUS_KM)
    kp, ku = pconstants.K_PENUMBRA, pconstants.K_UMBRA

    sin_angle_north = (sun_radius + kp) / shadow_dist
    sin_angle_south = (sun_radius - ku) / shadow_dist

    z_north = moon_z + (kp / sin_angle_north)
    z_south = moon_z - (ku / sin_angle_south)

    tangent_north = math.tan(math.asin(sin_angle_north))
    tangent_south = math.tan(math.asin(sin_angle_south))

    northern_limit = z_north * tangent_north
    southern_limit = z_south * tangent_south

    return (
        moon_x,
        moon_y,
        Angle(radians=shadow_decl).degrees,
        northern_limit,
        southern_limit,
        Angle(radians=sun_hour_angle).degrees,
        tangent_north,
        tangent_south,
    )


# ---------------------------------------------------------------------------
# First-Degree Polynomial Fit
# ---------------------------------------------------------------------------


def find_1st_degree_polynomial(
    val_0h: float,
    val_p1h: float,
    seconds: int = 3600,
) -> tuple[float, float, float, float]:
    """
    Compute a 1st-degree polynomial fit (linear) for a Besselian element.
    """
    H0 = val_0h
    H1 = val_p1h
    dH = (H1 - H0 + 180.0) % 360.0 - 180.0  # continuity
    mu = dH / (seconds / 3600.0)  # degrees per hour
    return H0, mu, 0.0, 0.0


# ---------------------------------------------------------------------------
# Third-Degree Polynomial Fit Using NumPy
# ---------------------------------------------------------------------------


def find_3rd_degree_polynomial(
    val_m2h: float,
    val_m1h: float,
    val_0h: float,
    val_p1h: float,
    val_p2h: float,
) -> tuple[float, float, float, float]:
    """
    Fit a cubic polynomial to five samples of a Besselian element using numpy.

    Returns coefficients (a0, a1, a2, a3) for a cubic polynomial.
    """
    # Matrix for times [-2, -1, 0, 1, 2]
    t = np.array([-2.0, -1.0, 0.0, 1.0, 2.0], dtype=float)
    b = np.array([val_m2h, val_m1h, val_0h, val_p1h, val_p2h], dtype=float)
    A = np.column_stack([np.ones_like(t), t, t * t, t * t * t])

    # Solve least-squares system
    coeffs, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

    return tuple(coeffs)
