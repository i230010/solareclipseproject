"""
pconstants.py
-----------------
Physical constants used throughout the solar eclipse computation modules.

Unless otherwise stated:
- All distances are expressed in kilometers.
- Earth-related constants follow the WGS-84 reference ellipsoid.
"""

# ---------------------------------------------------------------------------
# Solar constants
# ---------------------------------------------------------------------------

# Mean radius of the Sun in kilometers (IAU standard value)
SUN_RADIUS_KM: float = 695_700.0


# ---------------------------------------------------------------------------
# Earth constants (WGS-84)
# ---------------------------------------------------------------------------

# Equatorial radius of the Earth in kilometers
EARTH_RADIUS_KM: float = 6_378.137

# Square of Earth's first eccentricity (e^2), used in geodetic calculations
E_SQUARED: float = 0.006694385

# Complement of Earth's flattening factor: (1 - f)
# Used for converting between geodetic and geocentric coordinates
ONE_MINUS_F: float = 0.99664719

# Ellipsoid correction factor applied in geodetic-to-geographic conversions
ELLIPSOID_CORRECTION: float = 1.00336409


# ---------------------------------------------------------------------------
# Lunar constants
# ---------------------------------------------------------------------------

# Mean radius of the Moon in kilometers
MOON_RADIUS_KM: float = 1_737.400


# ---------------------------------------------------------------------------
# Eclipse-specific constants
# ---------------------------------------------------------------------------

# Longitude correction factor used in Delta-T-adjusted longitude calculations
# Units: degrees per second
DELTA_LAMBDA_FACTOR: float = 0.00417807

# Ratio of the Moon's apparent radius to the Sun's apparent radius
# for umbral eclipse calculations (Besselian elements)
K_UMBRA: float = 0.2722810

# Ratio of the Moon's apparent radius to the Sun's apparent radius
# for penumbral eclipse calculations (Besselian elements)
K_PENUMBRA: float = 0.2725076
