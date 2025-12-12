"""
psebessel.py
---------------

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


def km_to_earth_radii(km: float) -> float:
    """
    Convert a distance in kilometers to Earth radii.

    Parameters
    ----------
    km : float
        Distance in kilometers.

    Returns
    -------
    float
        Distance expressed in Earth radii.
    """
    return km / pconstants.EARTH_RADIUS_KM


def besselian_find(
    dt: pedatetime.datetime,
) -> tuple[float, float, float, float, float, float, float, float]:
    """
    Compute Besselian elements for a given datetime.

    These elements describe the geometry of a solar eclipse, including
    the position of the Moon's shadow relative to Earth and its limits.

    Parameters
    ----------
    dt : pedatetime.datetime
        The datetime for which to compute the Besselian elements.

    Returns
    -------
    tuple
        moon_x, moon_y : float
            Moon coordinates relative to shadow axis (Earth radii)
        shadow_decl_deg : float
            Shadow axis declination in degrees
        northern_limit, southern_limit : float
            Distances along shadow axis to northern and southern shadow limits
        sun_hour_angle_deg : float
            Hour angle of Sun minus shadow axis angle in degrees
        tangent_north, tangent_south : float
            Tangents for shadow limits
    """
    # Load planetary ephemeris and extract Earth, Sun, Moon
    planets = load(pdefilepath.EPHEM_PATH)
    earth, sun, moon = planets["earth"], planets["sun"], planets["moon"]

    # Convert datetime to Skyfield Time object in Terrestrial Time (TT)
    ts = load.timescale()
    ts.julian_calendar_cutoff = GREGORIAN_START
    t_sf = ts.tt(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
    
    # Compute apparent positions of Sun and Moon from Earth's center
    obs_sun = earth.at(t_sf).observe(sun)
    obs_moon = earth.at(t_sf).observe(moon)

    # Geocentric right ascension, declination, and distances
    sun_ra, sun_dec, sun_dist = obs_sun.radec()
    moon_ra, moon_dec, moon_dist = obs_moon.radec()

    # Convert distances to Earth radii
    sun_radius_earth_r = km_to_earth_radii(sun_dist.km)
    moon_radius_earth_r = km_to_earth_radii(moon_dist.km)

    # RA and Dec in radians
    sun_ra_rad = sun_ra.radians
    sun_dec_rad = sun_dec.radians
    moon_ra_rad = moon_ra.radians
    moon_dec_rad = moon_dec.radians

    # Rectangular coordinates of Sun and Moon in Earth radii
    sun_vec = vector.obj(
        x=sun_radius_earth_r * math.cos(sun_dec_rad) * math.cos(sun_ra_rad),
        y=sun_radius_earth_r * math.cos(sun_dec_rad) * math.sin(sun_ra_rad),
        z=sun_radius_earth_r * math.sin(sun_dec_rad),
    )
    moon_vec = vector.obj(
        x=moon_radius_earth_r * math.cos(moon_dec_rad) * math.cos(moon_ra_rad),
        y=moon_radius_earth_r * math.cos(moon_dec_rad) * math.sin(moon_ra_rad),
        z=moon_radius_earth_r * math.sin(moon_dec_rad),
    )

    # Vector from Moon to Sun and its magnitude
    shadow_vector = sun_vec - moon_vec
    shadow_distance = abs(shadow_vector)

    # Shadow axis angle in XY-plane
    shadow_axis_angle = math.atan2(shadow_vector.y, shadow_vector.x)

    # Shadow declination angle
    shadow_sin_decl = shadow_vector.z / shadow_distance
    shadow_decl = math.asin(shadow_sin_decl)

    # Hour angle of Sun minus shadow axis angle
    sun_hour_angle = (Angle(degrees=(t_sf.gast * 15)).radians - shadow_axis_angle) % (
        2 * math.pi
    )
    

    # Transform Moon coordinates relative to shadow axis
    moon_x = moon_radius_earth_r * (
        math.cos(moon_dec_rad) * math.sin(moon_ra_rad - shadow_axis_angle)
    )
    moon_y = moon_radius_earth_r * (
        (math.sin(moon_dec_rad) * math.cos(shadow_decl))
        - (
            math.cos(moon_dec_rad)
            * math.sin(shadow_decl)
            * math.cos(moon_ra_rad - shadow_axis_angle)
        )
    )
    moon_z = moon_radius_earth_r * (
        (math.sin(moon_dec_rad) * math.sin(shadow_decl))
        + (
            math.cos(moon_dec_rad)
            * math.cos(shadow_decl)
            * math.cos(moon_ra_rad - shadow_axis_angle)
        )
    )

    # Sun and Moon radii in Earth radii
    sun_radius = km_to_earth_radii(pconstants.SUN_RADIUS_KM)
    moon_radius = km_to_earth_radii(pconstants.MOON_RADIUS_KM)

    # Compute sine of angles for shadow limits
    sin_angle_north = (sun_radius + moon_radius) / shadow_distance
    sin_angle_south = (sun_radius - moon_radius) / shadow_distance

    # Z-coordinates for northern and southern shadow limits
    z_north = moon_z + (moon_radius / sin_angle_north)
    z_south = moon_z - (moon_radius / sin_angle_south)

    # Tangents for shadow limits
    tangent_north = math.tan(math.asin(sin_angle_north))
    tangent_south = math.tan(math.asin(sin_angle_south))

    # Distances along shadow axis to northern and southern limits
    northern_limit = z_north * tangent_north
    southern_limit = z_south * tangent_south

    return (
        moon_x,
        moon_y,
        Angle(radians=shadow_decl).degrees, #unfiltered micro
        northern_limit,
        southern_limit,
        Angle(radians=sun_hour_angle).degrees,
        tangent_north,
        tangent_south,
    )

def fit_hour_angle_polynomial(H_m2h, H_m1h, H_0h, H_p1h, H_p2h):
    """
    Fit a cubic polynomial H(t) = c0 + c1*t + c2*t^2 + c3*t^3
    to hour-angle values at t = -2h, -1h, 0h, +1h, +2h.

    Returns:
        coeffs : ndarray
            Polynomial coefficients [c0, c1, c2, c3]
        mu : float
            Linear coefficient c1 = μ (degrees per hour)
    """
    # Vandermonde matrix for t = -2, -1, 0, 1, 2
    A = np.array([
        [1, -2, 4, 8],
        [1, -1, 1, -1],
        [1,  0, 0, 0],
        [1,  1, 1, 1],
        [1,  2, 4, 8]
    ], dtype=float)

    # Hour-angle samples
    H_samples = np.array([H_m2h, H_m1h, H_0h, H_p1h, H_p2h], dtype=float)

    # Solve least-squares system
    coeffs, _, _, _ = np.linalg.lstsq(A, H_samples, rcond=None)

    # μ is the linear coefficient c1
    mu = coeffs[1]

    return coeffs, mu

def compute_mu(dt, dt_seconds=1):
    t0 = dt
    t1 = dt + pedatetime.timedelta(0, 0, 0, dt_seconds)

    # Compute hour angle difference for Sun relative to shadow axis
    _, _, _, _, _, H0, _, _ = besselian_find(t0)
    _, _, _, _, _, H1, _, _ = besselian_find(t1)

    # Ensure continuity across 360°
    dH = (H1 - H0 + 180) % 360 - 180

    # μ = change in degrees per hour
    mu = dH / (dt_seconds / 3600)  # degrees per hour
    return H0, mu


# -----------------------------------------------------------------------------
# Fit polynomial to Besselian elements
# -----------------------------------------------------------------------------
def find_besselian_polynomial(
    val_m2h: float, val_m1h: float, val_0h: float, val_p1h: float, val_p2h: float
) -> np.ndarray:
    """
    Fit a cubic polynomial to five samples of a Besselian element.

    Parameters
    ----------
    val_m2h, val_m1h, val_0h, val_p1h, val_p2h : float
        Besselian element values at times -2h, -1h, 0h, +1h, +2h.

    Returns
    -------
    coeffs : ndarray
        Array of 4 polynomial coefficients (cubic fit).
    """
    # Vandermonde-like matrix for cubic polynomial fit
    A = np.array(
        [
            [1, -2, 4, 8],
            [1, -1, 1, -1],
            [1, 0, 0, 0],
            [1, 1, 1, 1],
            [1, 2, 4, 8],
        ],
        dtype=float,
    )

    # Besselian element values at five sample times
    b = np.array([val_m2h, val_m1h, val_0h, val_p1h, val_p2h], dtype=float)

    # Solve least-squares system to get cubic polynomial coefficients
    coeffs, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

    return coeffs
