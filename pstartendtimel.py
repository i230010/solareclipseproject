"""
pstartendtimel.py
-----------------
Compute start and end times for penumbral or umbral intersections
using cubic Besselian polynomials and root-finding.
"""

from scipy.optimize import fsolve
from typing import Tuple, List


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


def penumbra_distance(
    t: float, x_coeffs: List[float], y_coeffs: List[float], l_coeffs: List[float]
) -> float:
    """
    Compute the penumbral distance function D(t):
        D(t) = sqrt(x(t)^2 + y(t)^2) - (1 + L(t))

    This function is used by fsolve to find when the penumbra or umbra
    first contacts or leaves the Earth's disk.

    Parameters
    ----------
    t : float
        Time parameter (TDT units).
    x_coeffs, y_coeffs, l_coeffs : list of float
        Cubic polynomial coefficients for X(t), Y(t), and L(t).

    Returns
    -------
    float
        Distance whose zero-crossings correspond to start/end times.
    """
    x = poly(x_coeffs, t)
    y = poly(y_coeffs, t)
    radius = poly(l_coeffs, t)
    return (x**2 + y**2) ** 0.5 - (1 + radius)


def startendtime(
    x_coeffs: List[float], y_coeffs: List[float], l_coeffs: List[float]
) -> Tuple[float, float]:
    """
    Solve for start and end times of the penumbra or umbra intersection.

    Uses fsolve with initial guesses to bracket the expected eclipse times.

    Parameters
    ----------
    x_coeffs, y_coeffs, l_coeffs : list of float
        Cubic polynomial coefficients for Besselian elements X(t), Y(t), L(t).

    Returns
    -------
    tuple of float
        (start_time, end_time) in the same units as the input polynomials.
    """
    # Find start time (initial guess -4 hours)
    te_start = fsolve(lambda t: penumbra_distance(t, x_coeffs, y_coeffs, l_coeffs), -4)
    # Find end time (initial guess +4 hours)
    te_end = fsolve(lambda t: penumbra_distance(t, x_coeffs, y_coeffs, l_coeffs), 4)

    # Safety: fsolve may return an array, take first element
    start_time = float(te_start[0]) if te_start.size > 0 else None
    end_time = float(te_end[0]) if te_end.size > 0 else None

    return start_time, end_time
