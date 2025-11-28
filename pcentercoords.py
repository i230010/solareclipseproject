"""
pcentercoords.py
-----------------
Computes geodetic latitude and longitude of the eclipse center
from Besselian elements using polynomial fits and Earth constants.
"""

import math
from skyfield.units import Angle
import pconstants


def poly(coeffs: list[float], t: float) -> float:
    """
    Evaluate a cubic polynomial at time t.

    P(t) = c0 + c1*t + c2*t^2 + c3*t^3

    Parameters
    ----------
    coeffs : list of float
        Polynomial coefficients [c0, c1, c2, c3]
    t : float
        Input value

    Returns
    -------
    float
        Polynomial evaluated at t
    """
    return coeffs[0] + coeffs[1] * t + coeffs[2] * t**2 + coeffs[3] * t**3


def coords(
    X_coeffs: list[float],
    Y_coeffs: list[float],
    declination_coeffs: list[float],
    longitude_coeffs: list[float],
    delta_t: float,
    t: float,
) -> tuple[float, float] | tuple[None, None]:
    """
    Compute geodetic latitude and longitude from Besselian elements.

    Parameters
    ----------
    X_coeffs : list of float
        Polynomial coefficients for X(t)
    Y_coeffs : list of float
        Polynomial coefficients for Y(t)
    declination_coeffs : list of float
        Polynomial coefficients for declination D(t) in degrees
    longitude_coeffs : list of float
        Polynomial coefficients for Moon's longitude Micro(t) in degrees
    delta_t : float
        Time difference correction (seconds)
    t : float
        Time parameter for polynomial evaluation

    Returns
    -------
    tuple
        (latitude, longitude) in degrees, or (None, None) if computation is invalid
    """

    # Evaluate polynomials
    X = poly(X_coeffs, t)
    Y = poly(Y_coeffs, t)
    decl_rad = Angle(degrees=poly(declination_coeffs, t)).radians
    longitude_rad = Angle(degrees=poly(longitude_coeffs, t)).radians

    # Correct for Earth's eccentricity
    omega = 1 / math.sqrt(1 - pconstants.E_SQUARED * (math.cos(decl_rad) ** 2))

    # Auxiliary variables for geodetic transformation
    Y_corr = omega * Y
    b1 = omega * math.sin(decl_rad)
    b2 = pconstants.ONE_MINUS_F * omega * math.cos(decl_rad)

    # Compute B for latitude calculation
    B_sq = 1 - X**2 - Y_corr**2
    if B_sq < 0:
        return None, None
    B = math.sqrt(B_sq)

    # Compute geodetic latitude
    sin_phi1 = B * b1 + Y_corr * b2
    phi1 = math.asin(sin_phi1)
    phi = math.atan(pconstants.ELLIPSOID_CORRECTION * math.tan(phi1))

    # Handle potential division by zero in hour angle calculation
    cos_phi1 = math.cos(phi1)
    if abs(cos_phi1) < 1e-12:
        return None, None

    sin_H = X / cos_phi1
    cos_H = (B * b2 - Y_corr * b1) / cos_phi1
    H = math.atan2(sin_H, cos_H)

    # Compute geodetic longitude with delta_t correction
    lambda_geo = (
        longitude_rad - H - (pconstants.DELTA_LAMBDA_FACTOR * delta_t * math.pi / 180)
    )

    # Convert to degrees and normalize to [-90, 90] and [-180, 180]
    lat_uncorrected = Angle(radians=phi).degrees
    lon_uncorrected = Angle(radians=lambda_geo).degrees * -1

    latitude = ((lat_uncorrected + 90) % 180) - 90
    longitude = ((lon_uncorrected + 180) % 360) - 180

    return latitude, longitude
