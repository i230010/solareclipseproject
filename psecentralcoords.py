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
from typing import List, Tuple

from skyfield.units import Angle

import pconstants


# ---------------------------------------------------------------------------
# Polynomial evaluation utility
# ---------------------------------------------------------------------------


def poly(coeffs: List[float], t: float) -> float:
    """
    Evaluate a cubic polynomial:
        P(t) = c0 + c1*t + c2*t^2 + c3*t^3

    Parameters
    ----------
    coeffs : list of float
        Polynomial coefficients [c0, c1, c2, c3].
    t : float
        Input variable (time).

    Returns
    -------
    float
        Polynomial evaluated at t.
    """
    return coeffs[0] + coeffs[1] * t + coeffs[2] * t * t + coeffs[3] * t * t * t


# ---------------------------------------------------------------------------
# Geodetic coordinate computation
# ---------------------------------------------------------------------------


def coords(
    x_coeffs: List[float],
    y_coeffs: List[float],
    d_coeffs: List[float],
    m_coeffs: List[float],
    delta_t: float,
    t: float,
) -> Tuple[float | None, float | None]:
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
        Time used for polynomial evaluation (time)

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
    X: float = poly(x_coeffs, t)
    Y: float = poly(y_coeffs, t)
    d_rad: float = Angle(degrees=poly(d_coeffs, t)).radians  # declination in radians
    m_rad: float = Angle(
        degrees=poly(m_coeffs, t)
    ).radians  # Greenwich hour angle in radians

    # -----------------------------------------------------------------------
    # Ellipsoid correction factors
    # -----------------------------------------------------------------------
    e2: float = pconstants.E_SQUARED  # Earth's eccentricity squared
    one_minus_f: float = pconstants.ONE_MINUS_F  # 1 - flattening factor

    omega: float = 1.0 / math.sqrt(1.0 - e2 * math.cos(d_rad) ** 2)

    y1: float = omega * Y
    b1: float = omega * math.sin(d_rad)
    b2: float = one_minus_f * omega * math.cos(d_rad)

    # -----------------------------------------------------------------------
    # Radial distance term (B)
    # -----------------------------------------------------------------------
    Bsq: float = 1.0 - X * X - y1 * y1
    if Bsq < 0.0:
        # Invalid geometry: eclipse shadow does not intersect Earth
        return None, None
    B: float = math.sqrt(Bsq)

    # -----------------------------------------------------------------------
    # Geocentric latitude (phi1)
    # -----------------------------------------------------------------------
    sin_phi1: float = B * b1 + y1 * b2
    phi1: float = math.asin(sin_phi1)

    # Convert geocentric latitude to geodetic latitude
    phi: float = math.atan(pconstants.ELLIPSOID_CORRECTION * math.tan(phi1))
    cos_phi1: float = math.cos(phi1)

    # -----------------------------------------------------------------------
    # Hour angle (H)
    # -----------------------------------------------------------------------
    sin_H: float = X / cos_phi1
    cos_H: float = (B * b2 - y1 * b1) / cos_phi1
    H: float = math.atan2(sin_H, cos_H)

    # -----------------------------------------------------------------------
    # Corrected longitude (lambda_geo)
    # -----------------------------------------------------------------------
    lambda_geo: float = (
        m_rad - H - pconstants.DELTA_LAMBDA_FACTOR * delta_t * math.pi / 180.0
    ) % (2.0 * math.pi)

    # Convert to degrees
    lat: float = Angle(radians=phi).degrees  # ty:ignore[invalid-assignment]
    lon: float = -Angle(radians=lambda_geo).degrees  # ty:ignore[unsupported-operator]

    # Normalize latitude to [-90, 90] and longitude to [-180, 180]
    lat: float = ((lat + 90.0) % 180.0) - 90.0
    lon: float = ((lon + 180.0) % 360.0) - 180.0

    return lat, lon
