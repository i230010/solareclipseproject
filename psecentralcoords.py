"""
psecentralcoords.py
------------------
Compute geodetic latitude and longitude of the eclipse shadow axis
using Besselian polynomial evaluations and Earth ellipsoid corrections.

The computation accounts for:
- Earth's oblateness
- Greenwich hour angle
- Ellipsoid-based latitude correction
"""

import math
from typing import Sequence, Tuple

from skyfield.units import Angle

import pconstants


# ---------------------------------------------------------------------------
# Polynomial evaluation utility
# ---------------------------------------------------------------------------

def poly(coeffs: Sequence[float], t: float) -> float:
    """
    Evaluate a cubic polynomial at time t.

        P(t) = a0 + a1*t + a2*t^2 + a3*t^3

    Parameters
    ----------
    coeffs : Sequence[float]
        Cubic polynomial coefficients [a0, a1, a2, a3].
    t : float
        Time parameter (e.g., Julian date or Besselian time).

    Returns
    -------
    float
        Polynomial value at time t.
    """
    a0, a1, a2, a3 = coeffs
    return a0 + a1 * t + a2 * t * t + a3 * t * t * t


# ---------------------------------------------------------------------------
# Geodetic coordinate computation
# ---------------------------------------------------------------------------

def coords(
    x_coeffs: Sequence[float],
    y_coeffs: Sequence[float],
    d_coeffs: Sequence[float],
    m_coeffs: Sequence[float],
    delta_t: float,
    t: float,
) -> Tuple[float, float]:
    """
    Compute geodetic latitude and longitude of the eclipse shadow axis.

    Parameters
    ----------
    x_coeffs, y_coeffs, d_coeffs, m_coeffs : Sequence[float]
        Cubic polynomial coefficients for Besselian elements X, Y, D
        (declination), and M (Greenwich hour angle).
    delta_t : float
        Delta-T correction in minutes used for longitude adjustment.
    t : float
        Time used for polynomial evaluation (e.g., Julian date or
        Besselian time).

    Returns
    -------
    Tuple[float, float]
        Geodetic latitude and longitude in degrees (lat, lon).

    Raises
    ------
    ValueError
        If the eclipse shadow does not intersect Earth.
    """
    # -----------------------------------------------------------------------
    # Evaluate Besselian polynomials at time t
    # -----------------------------------------------------------------------
    X = poly(x_coeffs, t)
    Y = poly(y_coeffs, t)
    d_rad = Angle(degrees=poly(d_coeffs, t)).radians  # declination in radians
    m_rad = Angle(degrees=poly(m_coeffs, t)).radians  # Greenwich hour angle in radians

    # -----------------------------------------------------------------------
    # Ellipsoid correction factors
    # -----------------------------------------------------------------------
    e2 = pconstants.E_SQUARED            # Earth's eccentricity squared
    one_minus_f = pconstants.ONE_MINUS_F  # 1 - flattening factor

    omega = 1.0 / math.sqrt(1.0 - e2 * math.cos(d_rad) ** 2)

    y1 = omega * Y
    b1 = omega * math.sin(d_rad)
    b2 = one_minus_f * omega * math.cos(d_rad)

    # -----------------------------------------------------------------------
    # Radial distance term (B)
    # -----------------------------------------------------------------------
    Bsq = 1.0 - X * X - y1 * y1
    if Bsq < 0.0:
        # Invalid geometry: eclipse shadow does not intersect Earth
        return None, None
    B = math.sqrt(Bsq)

    # -----------------------------------------------------------------------
    # Geocentric latitude (phi1)
    # -----------------------------------------------------------------------
    sin_phi1 = B * b1 + y1 * b2
    phi1 = math.asin(sin_phi1)

    # Convert geocentric latitude to geodetic latitude
    phi = math.atan(pconstants.ELLIPSOID_CORRECTION * math.tan(phi1))
    cos_phi1 = math.cos(phi1)

    # -----------------------------------------------------------------------
    # Hour angle (H)
    # -----------------------------------------------------------------------
    sin_H = X / cos_phi1
    cos_H = (B * b2 - y1 * b1) / cos_phi1
    H = math.atan2(sin_H, cos_H)

    # -----------------------------------------------------------------------
    # Corrected longitude (lambda_geo)
    # -----------------------------------------------------------------------
    lambda_geo = (
        m_rad
        - H
        - pconstants.DELTA_LAMBDA_FACTOR * delta_t * math.pi / 180.0
    ) % (2.0 * math.pi)

    # Convert to degrees
    lat = Angle(radians=phi).degrees
    lon = -Angle(radians=lambda_geo).degrees

    # Normalize latitude to [-90, 90] and longitude to [-180, 180]
    lat = ((lat + 90.0) % 180.0) - 90.0
    lon = ((lon + 180.0) % 360.0) - 180.0

    return lat, lon
