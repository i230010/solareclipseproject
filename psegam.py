"""
psegam.py
-----------------
Computes the eclipse gamma parameter for a given moment in time.
"""

import math


def poly(coeffs: list, t: float) -> float:
    """
    Evaluate a cubic polynomial at time t.

    Polynomial form: a0 + a1*t + a2*t^2 + a3*t^3

    Args:
        coeffs (list[float]): Coefficients [a0, a1, a2, a3]
        t (float): Independent variable (time)

    Returns:
        float: Polynomial value at t
    """
    return coeffs[0] + coeffs[1] * t + coeffs[2] * t**2 + coeffs[3] * t**3


def gamma(x_coeffs: list, y_coeffs: list, t_max: float) -> float:
    """
    Compute the eclipse gamma parameter at a given time.

    Gamma represents the distance (in Earth radii) of the Moon's shadow
    axis from Earth's center. Sign convention: negative if the shadow axis
    is south of the equator (Y < 0).

    Args:
        x_coeffs (list[float]): Cubic coefficients for X Besselian element
        y_coeffs (list[float]): Cubic coefficients for Y Besselian element
        t_max (float): Time in Besselian time units

    Returns:
        float: Gamma parameter (signed distance from central eclipse line)
    """
    # Evaluate X and Y Besselian coordinates at t_max
    x_val = poly(x_coeffs, t_max)
    y_val = poly(y_coeffs, t_max)

    # Compute Euclidean distance from central eclipse line using math.hypot
    gamma_val = math.hypot(x_val, y_val)

    # Apply sign convention: negative if Y coordinate is south
    if y_val < 0:
        gamma_val = -gamma_val

    return gamma_val
