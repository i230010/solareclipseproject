"""
pconstants.py
-----------------

This module stores planet-related physical constants used throughout
the project. Values are taken from standard astronomical and geodetic
references.

"""

# Radius of the Sun in kilometers
SUN_RADIUS_KM: float = 695_700.0

# Equatorial radius of the Earth in kilometers
EARTH_RADIUS_KM: float = 6_378.137

# Mean radius of the Moon in kilometers
MOON_RADIUS_KM: float = 1_737.400

# Earth's eccentricity squared (WGS-84)
E_SQUARED: float = 0.006694385

# Ellipsoid correction factor used for geodetic calculations
ELLIPSOID_CORRECTION: float = 1.00336409
