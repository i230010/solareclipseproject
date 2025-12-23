"""
psestartendtime.py
-----------------
Compute start and end times for penumbral or umbral intersections
using cubic Besselian polynomials and numerical root-finding.

The zero crossings of the distance function correspond to the
first and last contact of the eclipse shadow with the Earth.
"""

from typing import Sequence, Tuple
import math

from scipy.optimize import brentq


# ---------------------------------------------------------------------------
# Polynomial utilities
# ---------------------------------------------------------------------------


def poly(coeffs: Sequence[float], t: float) -> float:
    """
    Evaluate a cubic polynomial:

        P(t) = c0 + c1*t + c2*t^2 + c3*t^3

    Parameters
    ----------
    coeffs : sequence of float
        Cubic polynomial coefficients [c0, c1, c2, c3].
    t : float
        Time variable.

    Returns
    -------
    float
        Polynomial value at time t.
    """
    c0, c1, c2, c3 = coeffs
    return c0 + c1 * t + c2 * t * t + c3 * t * t * t


# ---------------------------------------------------------------------------
# Eclipse distance function
# ---------------------------------------------------------------------------


def penumbra_distance(
    t: float,
    x_coeffs: Sequence[float],
    y_coeffs: Sequence[float],
    l_coeffs: Sequence[float],
) -> float:
    """
    Compute the distance function whose zero defines shadow contact.

        D(t) = sqrt(x(t)^2 + y(t)^2) - (1 + L(t))

    A root of D(t) corresponds to the time when the penumbral or umbral
    boundary intersects the Earth's disk.

    Parameters
    ----------
    t : float
        Time parameter (same units as the Besselian polynomials).
    x_coeffs, y_coeffs, l_coeffs : sequence of float
        Cubic polynomial coefficients for X(t), Y(t), and L(t).

    Returns
    -------
    float
        Signed distance value.
    """
    x = poly(x_coeffs, t)
    y = poly(y_coeffs, t)
    radius = poly(l_coeffs, t)

    # hypot(x, y) computes sqrt(x^2 + y^2) in a numerically stable way
    return math.hypot(x, y) - (1.0 + radius)


# ---------------------------------------------------------------------------
# Start and end time solver
# ---------------------------------------------------------------------------


def startendtime(
    x_coeffs: Sequence[float],
    y_coeffs: Sequence[float],
    l_coeffs: Sequence[float],
    t_start: float = -6.0,
    t_mid: float = 0.0,
    t_end: float = 6.0,
) -> Tuple[float, float]:
    """
    Solve for the start and end times of penumbral or umbral contact.

    This function assumes the eclipse event is centered near t = 0 and
    that the distance function changes sign across the provided brackets.

    Parameters
    ----------
    x_coeffs, y_coeffs, l_coeffs : sequence of float
        Cubic polynomial coefficients for Besselian elements X(t), Y(t), L(t).
    t_start : float, optional
        Lower bound for the start time search (default: -6).
    t_mid : float, optional
        Midpoint separating start and end roots (default: 0).
    t_end : float, optional
        Upper bound for the end time search (default: 6).

    Returns
    -------
    tuple of float
        (start_time, end_time) in the same units as the input polynomials.

    Raises
    ------
    ValueError
        If the root is not bracketed within the provided intervals.
    """
    # Solve for first contact (ingress)
    start_time = brentq(
        penumbra_distance,
        t_start,
        t_mid,
        args=(x_coeffs, y_coeffs, l_coeffs),
    )

    # Solve for last contact (egress)
    end_time = brentq(
        penumbra_distance,
        t_mid,
        t_end,
        args=(x_coeffs, y_coeffs, l_coeffs),
    )

    return start_time, end_time
