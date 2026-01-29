# ppartialmazcoords.py
# Module for calculating the magnitude of a partial solar eclipse
# and for finding the geographic coordinates of maximum partial eclipse
# using Besselian elements and numerical optimization.
from typing import Tuple, List

import pselocalcirumstances  # Custom module to compute local eclipse circumstances
from scipy.optimize import minimize  # For numerical optimization


def magnitude(
    Xa: float,
    Ya: float,
    Da: float,
    Ma: float,
    L1a: float,
    L2a: float,
    tanf1: float,
    tanf2: float,
    T0: float,
    lat: float,
    lon: float,
    height: float,
    delta_t: float,
) -> float:
    """
    Compute the magnitude of a partial solar eclipse at a specific observer location.

    Parameters:
        Xa, Ya, Da, Ma, L1a, L2a, tanf1, tanf2, T0, delta_t: Besselian elements
        lat, lon: observer latitude and longitude in degrees
        height: observer height in metre

    Returns:
        magnitude (float): The fraction of the Sun's diameter obscured by the Moon
    """

    # Use the "unsafe" local circumstances function from plocalcirumstances module.
    # Note: longitude is negated here to match the expected convention in that module.
    e = pselocalcirumstances.get_local_circumstances_unsafe(
        Xa, Ya, Da, Ma, L1a, L2a, tanf1, tanf2, T0, lat, -lon, height, delta_t
    )

    # Return the magnitude as a float
    return float(e["magnitude"])


def coords(
    Xa: float,
    Ya: float,
    Da: float,
    Ma: float,
    L1a: float,
    L2a: float,
    tanf1: float,
    tanf2: float,
    T0: float,
    delta_t: float,
    height: float = 0.0,
) -> Tuple[float, float]:
    """
    Find the latitude and longitude of maximum partial solar eclipse
    using numerical optimization (scipy.optimize.minimize).

    Parameters:
        Xa, Ya, Da, Ma, L1a, L2a, tanf1, tanf2, T0, delta_t: Besselian elements
        height: observer height in km (default 0)

    Returns:
        best_lat (float): Latitude of maximum partial eclipse
        best_lon (float): Longitude of maximum partial eclipse
        max_magnitude (float): Magnitude at the optimal location
    """

    # Define a helper function to optimize.
    # We minimize the negative of magnitude to perform maximization.
    def mag_to_optimize(coords_array: Tuple[float, float]) -> float:
        lat, lon = coords_array  # Unpack latitude and longitude from optimization array
        # Return negative magnitude for minimization
        return -magnitude(
            Xa, Ya, Da, Ma, L1a, L2a, tanf1, tanf2, T0, lat, lon, height, delta_t
        )

    # Set bounds for latitude (-90° to 90°) and longitude (-180° to 180°)
    bounds: List[Tuple[float, float]] = [(-90, 90), (-180, 180)]

    # Initial guess for optimization (equatorial point)
    x0: List[float] = [0, 0]

    # Perform optimization using L-BFGS-B method, suitable for bounded problems
    # ftol controls the tolerance for convergence
    result = minimize(
        mag_to_optimize, x0, bounds=bounds, method="L-BFGS-B", options={"ftol": 1e-10}
    )

    # Extract optimized latitude and longitude
    best_lat, best_lon = result.x

    # Return latitude, longitude
    return best_lat, best_lon
