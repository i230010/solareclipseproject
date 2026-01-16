"""
psegam.py
---------
Compute the eclipse gamma parameter at a given moment in time.

The gamma parameter represents the signed distance (in Earth radii)
between the Moon's shadow axis and the Earth's center. By convention,
gamma is negative when the shadow axis lies south of the Earth's equator.
"""

import math
from typing import List


# ---------------------------------------------------------------------------
# Polynomial utilities
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
# Gamma computation
# ---------------------------------------------------------------------------


def gamma(
    x_coeffs: List[float],
    y_coeffs: List[float],
    t_max: float,
) -> float:
    """
    Compute the eclipse gamma parameter.

    Gamma is the signed distance (in Earth radii) of the Moon's shadow axis
    from the Earth's center at the specified time.

    Sign convention:
    - gamma < 0 -> shadow axis south of the equator (Y < 0)
    - gamma > 0 -> shadow axis north of the equator (Y >= 0)

    Parameters
    ----------
    x_coeffs : sequence of float
        Cubic coefficients for the X Besselian element.
    y_coeffs : sequence of float
        Cubic coefficients for the Y Besselian element.
    t_max : float
        Time of maximum eclipse in Besselian time units.

    Returns
    -------
    float
        Signed gamma parameter.
    """
    # Evaluate X and Y Besselian elements at t_max
    x_val: float = poly(x_coeffs, t_max)
    y_val: float = poly(y_coeffs, t_max)

    # Compute distance from Earth's center using a stable Euclidean norm
    gamma_val: float = math.hypot(x_val, y_val)

    # Apply sign convention based on the Y coordinate
    return -gamma_val if y_val < 0.0 else gamma_val
