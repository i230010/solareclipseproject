"""
_besselians.py
---------------

Computes Besselian elements for a given datetime and fits
polynomials to these elements for smooth interpolation.

"""

import math
import numpy as np
from skyfield.api import load
from skyfield.units import Angle
import vector

import pconstants
import pdefilepath
import pedatetime


def kmtoearthradii(x: float) -> float:
    """
    Convert a distance in kilometers to Earth radii.

    Parameters
    ----------
    x : float
        Distance in kilometers.

    Returns
    -------
    float
        Distance expressed in Earth radii.
    """
    return x / pconstants.EARTH_RADIUS_KM


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
        x, y : float
            Moon coordinates relative to shadow axis (Earth radii)
        d : float
            Shadow axis declination in degrees
        l1, l2 : float
            Distances along shadow axis to northern and southern shadow limits
        u : float
            Hour angle of Sun minus shadow axis angle in degrees
        tanf1, tanf2 : float
            Tangents for shadow limits
    """
    # Load planetary ephemeris and extract Earth, Sun, Moon
    planets = load(pdefilepath.EPHEM_PATH)
    earth, sun, moon = planets["earth"], planets["sun"], planets["moon"]

    # Convert datetime to Skyfield Time object
    ts = load.timescale()
    t = ts.utc(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)

    # Compute apparent positions of Sun and Moon as seen from Earth
    obs_sun = earth.at(t).observe(sun)
    obs_moon = earth.at(t).observe(moon)

    # Geocentric right ascension, declination, and distances
    sra, sdec, sdist = obs_sun.radec()
    mra, mdec, mdist = obs_moon.radec()

    # Convert distances to Earth radii
    rp = kmtoearthradii(sdist.km)
    r = kmtoearthradii(mdist.km)

    # RA and Dec in radians
    alp = sra.radians
    dlp = sdec.radians
    al = mra.radians
    dl = mdec.radians

    # Rectangular coordinates of Sun and Moon in Earth radii
    xs = rp * math.cos(dlp) * math.cos(alp)
    ys = rp * math.cos(dlp) * math.sin(alp)
    zs = rp * math.sin(dlp)

    xm = r * math.cos(dl) * math.cos(al)
    ym = r * math.cos(dl) * math.sin(al)
    zm = r * math.sin(dl)

    # Create vector objects for Sun and Moon
    Sar = vector.obj(x=xs, y=ys, z=zs)
    Mar = vector.obj(x=xm, y=ym, z=zm)

    # Vector from Moon to Sun and its magnitude
    Gar = Sar - Mar
    G = abs(Gar)

    # Shadow axis angle in XY-plane
    tana = (Gar.y / G) / (Gar.x / G)
    a = math.atan(tana)

    # Sine of shadow declination
    sind = Gar.z / G

    # Hour angle of Sun minus shadow axis angle
    u = Angle(degrees=(t.gast * 15)).radians - a

    # Shadow axis declination
    d = math.asin(sind)

    # Transform Moon coordinates relative to shadow axis
    x = r * (math.cos(dl) * math.sin(al - a))
    y = r * (
        (math.sin(dl) * math.cos(d)) - (math.cos(dl) * math.sin(d) * math.cos(al - a))
    )
    z = r * (
        (math.sin(dl) * math.sin(d)) + (math.cos(dl) * math.cos(d) * math.cos(al - a))
    )

    # Sun and Moon radii in Earth radii
    ds = kmtoearthradii(pconstants.SUN_RADIUS_KM)
    k = kmtoearthradii(pconstants.MOON_RADIUS_KM)

    # Compute sine of angles for shadow limits
    sinf1 = (ds + k) / G
    sinf2 = (ds - k) / G

    # Z-coordinates for northern and southern shadow limits
    c1 = z + (k / sinf1)
    c2 = z - (k / sinf2)

    # Tangents for shadow limits
    tanf1 = math.tan(math.asin(sinf1))
    tanf2 = math.tan(math.asin(sinf2))

    # Distances along shadow axis to northern and southern limits
    l1 = c1 * tanf1
    l2 = c2 * tanf2

    return (
        x,
        y,
        Angle(radians=d).degrees,
        l1,
        l2,
        Angle(radians=u).degrees,
        tanf1,
        tanf2,
    )


# -----------------------------------------------------------------------------
# Fit polynomial to Besselian elements
# -----------------------------------------------------------------------------
def find_besselian_polynomial(
    tm2: tuple, tm1: tuple, t: tuple, tp1: tuple, tp2: tuple
) -> tuple:
    """
    Fit a cubic polynomial to five samples of a Besselian element.

    Parameters
    ----------
    tm2, tm1, t, tp1, tp2 : float
        Values at times -2h, -1h, 0h, +1h, +2h.

    Returns
    -------
    coeffs : ndarray
        Array of 4 polynomial coefficients (cubic fit).
    """
    # Vandermonde-like matrix for polynomial fit of degree 3
    A = np.array(
        [[1, -2, 4, 8], [1, -1, 1, -1], [1, 0, 0, 0], [1, 1, 1, 1], [1, 2, 4, 8]],
        dtype=float,
    )

    # Besselian element values at five sample times
    b = np.array([tm2, tm1, t, tp1, tp2], dtype=float)

    # Solve least-squares system to get polynomial coefficients
    coeffs, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

    return coeffs
