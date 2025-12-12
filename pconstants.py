"""
pconstants.py
-----------------
Stores physical constants used throughout the eclipse computations.

All distances are in kilometers unless specified otherwise.
"""

# Radius of the Sun in kilometers
SUN_RADIUS_KM: float = 695_700.0

# Equatorial radius of the Earth in kilometers
EARTH_RADIUS_KM: float = 6_378.137

# Mean radius of the Moon in kilometers
MOON_RADIUS_KM: float = 1_737.400

# Earth's eccentricity squared (WGS-84) used for geodetic calculations
E_SQUARED: float = 0.006694385

# Ellipsoid correction factor for converting between geodetic and geographic coordinates
ELLIPSOID_CORRECTION: float = 1.00336409

# Complement of Earth's flattening factor (1 - f)
ONE_MINUS_F: float = 0.99664719

# Correction factor used in longitude calculation for delta_t adjustment
DELTA_LAMBDA_FACTOR: float = 0.00417807  # degrees per second

CENTRAL_GAMMA_LIMIT: float = 0.9972

EXIST_GAMMA_LIMIT: float = 1.56
