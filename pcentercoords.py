"""
pcentercoords.py
-----------------
Computes geodetic latitude and longitude using polynomial evaluations and Earth ellipsoid corrections.

Uses polynomial coefficients and time to calculate the coordinates, adjusting for Greenwich hour angle
and Earth's oblateness.
"""

import math
from skyfield.units import Angle
import pconstants


def poly(a: list, t: float) -> float:
    """
    Evaluates a cubic polynomial at a specific time.

    Parameters:
        a (list): Polynomial coefficients [a0, a1, a2, a3] for the polynomial equation.
        t (float): The time at which to evaluate the polynomial. Typically the Julian Date or similar time.

    Returns:
        float: The value of the polynomial at time t.
    """
    # Polynomial evaluation: a0 + a1*t + a2*t^2 + a3*t^3
    return a[0] + a[1] * t + a[2] * t**2 + a[3] * t**3


def coords(Xa: list, Ya: list, Da: list, Ma: list, delta_t: float, T: float) -> tuple:
    """
    Computes geodetic latitude and longitude using polynomial coefficients for
    the astronomical components and the time correction.

    This function applies the Besselian evaluation and Earth ellipsoid corrections
    to compute the geodetic coordinates, adjusting for the oblateness of Earth
    and correcting for the Greenwich hour angle.

    Parameters:
        Xa, Ya, Da, Ma (list): Polynomial coefficients for X, Y, D (Declination), and M (Greenwich hour angle).
        delta_t (float): Time correction factor (in minutes) used for adjusting the longitude.
        T (float): Time used for polynomial evaluations (e.g., Julian Date).

    Returns:
        tuple: Geodetic latitude and longitude in degrees.
               If an error occurs in calculation (e.g., negative Bsq), returns (None, None).
    """
    # Perform polynomial evaluation for X, Y, D, and M at time T
    X = poly(Xa, T)
    Y = poly(Ya, T)
    d = Angle(degrees=poly(Da, T)).radians  # Declination in radians
    m = Angle(degrees=poly(Ma, T)).radians  # Greenwich hour angle in radians

    # Earth ellipsoid constants
    e2 = pconstants.E_SQUARED  # Square of Earth's eccentricity
    one_minus_f = pconstants.ONE_MINUS_F  # Earth flattening parameter

    # Ellipsoid corrections
    # Omega is the inverse of the square root of (1 - e2 * cos(d)^2)
    omega = 1.0 / math.sqrt(1 - e2 * (math.cos(d) ** 2))
    y1 = omega * Y  # Y' component, scaled by omega
    b1 = omega * math.sin(d)  # Sine of declination, scaled by omega
    b2 = (
        one_minus_f * omega * math.cos(d)
    )  # Cosine of declination, scaled by one_minus_f * omega

    # Compute B, the radial distance correction term
    Bsq = 1 - X**2 - y1**2  # Bsq = 1 - X^2 - Y'^2
    if Bsq < 0:
        return None, None
    B = math.sqrt(Bsq)  # B = sqrt(1 - X^2 - Y'^2)

    # Latitude components
    # Calculate the sin of the geocentric latitude (phi1) using the components of B
    sinphi1 = B * b1 + y1 * b2
    phi1 = math.asin(sinphi1)  # Geocentric latitude in radians

    # Oblateness correction to convert from geocentric to geodetic latitude
    # Geodetic latitude is corrected by Earth's oblateness factor (ELLIPSOID_CORRECTION)
    phi = math.atan(
        pconstants.ELLIPSOID_CORRECTION * math.tan(phi1)
    )  # Geodetic latitude in radians

    # Hour angle (H)
    sinH = X / math.cos(phi1)  # Sine of the hour angle
    cosH = (B * b2 - y1 * b1) / math.cos(phi1)  # Cosine of the hour angle
    H = math.atan2(sinH, cosH)  # Hour angle in radians

    # Compute the corrected longitude (lambda_geo)
    lambda_geo = (m - H - pconstants.DELTA_LAMBDA_FACTOR * delta_t * math.pi / 180) % (
        2 * math.pi
    )

    # Convert geodetic latitude and longitude from radians to degrees
    lat_uncorrected, lon_uncorrected = (
        Angle(radians=phi).degrees,
        (Angle(radians=lambda_geo).degrees * -1),
    )

    # Longitude is adjusted to fall within the range [-180, 180]
    lat = ((lat_uncorrected + 90) % 180) - 90
    lon = ((lon_uncorrected + 180) % 360) - 180
    return lat, lon
