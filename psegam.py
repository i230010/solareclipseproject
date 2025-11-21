"""
psegam.py
-----------------

Computes the eclipse gamma parameter for a given moment in time.

"""

from skyfield.api import load
import vector
import math

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


def gamma(dt: pedatetime.datetime) -> float:
    """
    Compute the gamma parameter for a solar eclipse at a given datetime.

    Gamma is the distance of the Moon's shadow axis from Earth's center
    (in Earth radii) at the instant of greatest eclipse. Positive gamma
    means the shadow passes north of Earth's center, negative means south.

    Parameters
    ----------
    dt : pedatetime.datetime
        The datetime at which to compute the gamma value.

    Returns
    -------
    float
        The gamma value in Earth radii.
    """
    # Load the timescale for astronomical calculations
    ts = load.timescale()
    # Convert the input datetime to a Skyfield Time object
    t = ts.utc(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)

    # Load planetary ephemeris from the DE file
    planets = load(pdefilepath.EPHEM_PATH)
    earth, sun, moon = planets["earth"], planets["sun"], planets["moon"]

    # Compute apparent positions of Sun and Moon as seen from Earth
    obs_sun = earth.at(t).observe(sun)
    obs_moon = earth.at(t).observe(moon)

    # Get right ascension (RA), declination (Dec), and distance for Sun and Moon
    sra, sdec, sdist = obs_sun.radec()
    mra, mdec, mdist = obs_moon.radec()

    # Convert distances from km to Earth radii
    rp = kmtoearthradii(sdist.km)  # Sun distance
    r = kmtoearthradii(mdist.km)  # Moon distance

    # Extract RA and Dec in radians
    alp = sra.radians
    dlp = sdec.radians
    al = mra.radians
    dl = mdec.radians

    # Convert spherical coordinates (RA, Dec, distance) to Cartesian coordinates
    xs = rp * math.cos(dlp) * math.cos(alp)
    ys = rp * math.cos(dlp) * math.sin(alp)
    zs = rp * math.sin(dlp)

    xm = r * math.cos(dl) * math.cos(al)
    ym = r * math.cos(dl) * math.sin(al)
    zm = r * math.sin(dl)

    # Create vector objects for Sun and Moon positions
    Sar = vector.obj(x=xs, y=ys, z=zs)
    Mar = vector.obj(x=xm, y=ym, z=zm)

    # Compute the vector from Moon to Sun
    Gar = Sar - Mar
    # Unit vector along Moon-to-Sun direction (shadow axis)
    G = Gar.unit()

    # Compute the projection of Moon vector onto shadow axis
    u = -(Mar.dot(G))
    p = Mar + G * u

    # Magnitude of the perpendicular distance from Earth's center to shadow axis
    D = p.mag
    # Determine sign based on whether the shadow passes north (+) or south (-) of Earth's center
    sign = +1 if p.z > 0 else -1
    gamma = sign * D

    return gamma
