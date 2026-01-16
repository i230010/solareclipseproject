"""
plocalcircumstances.py

Local solar eclipse circumstance calculations based on Besselian elements.

This module computes fundamental eclipse arguments, local contact times,
sun/moon positions, and eclipse magnitudes for a given observer location.

NOTE:
This code is borrowed from the gmiller123456 Solar Eclipse Viewer project.
It has been commented for clarity, but its functionality has not been changed.
"""

import math
from typing import List
from skyfield.units import Angle
import pconstants

def poly(coeffs: List[float], t: float) -> float:
    """
    Evaluate a cubic polynomial:

        P(t) = c0 + c1*t + c2*t^2 + c3*t^3

    Used for evaluating Besselian element series.
    """
    return coeffs[0] + coeffs[1] * t + coeffs[2] * t * t + coeffs[3] * t * t * t


def polyprime(coeffs: List[float], t: float) -> float:
    """
    First derivative of a cubic polynomial with respect to t.
    """
    return coeffs[1] + 2 * coeffs[2] * t + 3 * coeffs[3] * t * t

def to_m(km: float) -> float:
    return km*1000

def get_fundamental_arguments_unsafe(
    Xa, Ya, Da, Ma, L1a, L2a, tanf1, tanf2, t, lat, lon, height, delta_t
) -> dict:
    """
    Compute fundamental Besselian arguments for a given observer and time.

    This version assumes longitude sign is already handled correctly.
    """
    rad: float = Angle(degrees=1).radians

    # Dictionary holding all computed parameters
    o: dict = {}

    # Moon shadow coordinates and declination
    o["X"] = poly(Xa, t)
    o["Y"] = poly(Ya, t)
    o["d"] = poly(Da, t)
    d1: float = Da[1]

    # Hour angle of the Sun
    o["M"] = poly(Ma, t)
    M1: float = Ma[1]

    # Time derivatives
    o["Xp"] = polyprime(Xa, t)
    o["Yp"] = polyprime(Ya, t)

    # Penumbral and umbral radii
    o["L1"] = poly(L1a, t)
    o["L2"] = poly(L2a, t)

    # Local hour angle
    o["H"] = o["M"] - lon - pconstants.DELTA_LAMBDA_FACTOR * delta_t

    # Geocentric latitude correction
    u1: float = math.atan(pconstants.ONE_MINUS_F * math.tan(lat * rad)) / rad
    rho_sin: float = pconstants.ONE_MINUS_F * math.sin(u1 * rad) + height / to_m(pconstants.EARTH_RADIUS_KM) * math.sin(lat * rad)
    rho_cos: float = math.cos(u1 * rad) + height / to_m(pconstants.EARTH_RADIUS_KM) * math.cos(lat * rad)

    # Observer coordinates projected onto fundamental plane
    o["xi"]: float = rho_cos * math.sin(o["H"] * rad)
    o["eta"]: float = rho_sin * math.cos(o["d"] * rad) - rho_cos * math.cos(
        o["H"] * rad
    ) * math.sin(o["d"] * rad)
    o["zeta"]: float = rho_sin * math.sin(o["d"] * rad) + rho_cos * math.cos(
        o["H"] * rad
    ) * math.cos(o["d"] * rad)

    # Time derivatives of observer coordinates
    o["xip"]: float = rad * M1 * rho_cos * math.cos(o["H"] * rad)
    o["etap"]: float = rad * (M1 * o["xi"] * math.sin(o["d"] * rad) - o["zeta"] * d1)

    # Corrected shadow radii
    o["L1p"]: float = o["L1"] - o["zeta"] * tanf1
    o["L2p"]: float = o["L2"] - o["zeta"] * tanf2

    # Relative shadow position and velocity
    o["u"]: float = o["X"] - o["xi"]
    o["v"]: float = o["Y"] - o["eta"]
    o["a"]: float = o["Xp"] - o["xip"]
    o["b"]: float = o["Yp"] - o["etap"]

    # Shadow velocity magnitude
    o["n"]: float = math.sqrt(o["a"] ** 2 + o["b"] ** 2)

    return o


def get_fundamental_arguments(
    Xa, Ya, Da, Ma, L1a, L2a, tanf1, tanf2, t, lat, lon, height, delta_t
) -> dict:
    """
    Safe wrapper that flips longitude sign convention.
    """
    return get_fundamental_arguments_unsafe(
        Xa, Ya, Da, Ma, L1a, L2a, tanf1, tanf2, t, lat, -lon, height, delta_t
    )


def sun_moon_pos(Xa, Ya, Da, Ma, L1a, L2a, tanf1, tanf2, t, lat, lon, height, delta_t):
    """
    Compute apparent Sun and Moon altitude/azimuth for a given time.

    Returns:
        sun_alt, sun_az, moon_alt, moon_az (degrees)

    Azimuth convention:
        0° = North, 90° = East
    """
    # Compute fundamental eclipse arguments
    o = get_fundamental_arguments_unsafe(
        Xa, Ya, Da, Ma, L1a, L2a, tanf1, tanf2, t, lat, -lon, height, delta_t
    )

    # Convert angles to radians
    lat_r: float = Angle(degrees=lat).radians
    d: float = Angle(degrees=o["d"]).radians
    H: float = Angle(degrees=o["H"]).radians

    # -------------------------
    # Sun altitude and azimuth
    # -------------------------
    sin_hs: float = math.sin(lat_r) * math.sin(d) + math.cos(lat_r) * math.cos(d) * math.cos(H)
    sun_alt: float = math.asin(sin_hs)

    sun_az: float = math.atan2(
        math.sin(H), math.cos(H) * math.sin(lat_r) - math.tan(d) * math.cos(lat_r)
    )
    sun_az: float = (sun_az + math.pi) % (2 * math.pi)

    # -------------------------
    # Moon position relative to Sun
    # -------------------------
    u: float = o["u"]
    v: float = o["v"]

    moon_dec: float = d + v
    moon_H: float = H - u / math.cos(d)

    # -------------------------
    # Moon altitude and azimuth
    # -------------------------
    sin_hm: float = math.sin(lat_r) * math.sin(moon_dec) + math.cos(lat_r) * math.cos(
        moon_dec
    ) * math.cos(moon_H)
    moon_alt: float = math.asin(sin_hm)

    moon_az: float = math.atan2(
        math.sin(moon_H),
        math.cos(moon_H) * math.sin(lat_r) - math.tan(moon_dec) * math.cos(lat_r),
    )
    moon_az: float = (moon_az + math.pi) % (2 * math.pi)

    return (
        Angle(radians=sun_alt).degrees,
        Angle(radians=sun_az).degrees,
        Angle(radians=moon_alt).degrees,
        Angle(radians=moon_az).degrees,
    )


def get_local_circumstances_unsafe(
    Xa, Ya, Da, Ma, L1a, L2a, tanf1, tanf2, T0, lat, lon, height, delta_t
):
    """
    Compute local eclipse circumstances using iterative minimization.
    """
    rad: float = Angle(degrees=1).radians
    t: float = 0.0
    taum: float = 1e9
    iterations: int = 0
    MAX_ITERATIONS: int = 50

    # Iteratively solve for time of closest approach
    while abs(taum) > 1e-5 and iterations < MAX_ITERATIONS:
        o = get_fundamental_arguments_unsafe(
            Xa, Ya, Da, Ma, L1a, L2a, tanf1, tanf2, t, lat, lon, height, delta_t
        )
        taum: float = -(o["u"] * o["a"] + o["v"] * o["b"]) / (o["n"] ** 2)
        t += taum
        iterations += 1

    # Eclipse magnitude and obscuration ratio
    m: float = math.sqrt(o["u"] ** 2 + o["v"] ** 2)
    G: float = (o["L1p"] - m) / (o["L1p"] + o["L2p"])
    A: float = (o["L1p"] - o["L2p"]) / (o["L1p"] + o["L2p"])

    # Solar altitude and azimuth at maximum eclipse
    sinh: float = math.sin(o["d"] * rad) * math.sin(lat * rad) + math.cos(
        o["d"] * rad
    ) * math.cos(lat * rad) * math.cos(o["H"] * rad)
    h: float = math.asin(sinh) / rad
    q: float = math.asin(math.cos(lat * rad) * math.sin(o["H"] * rad) / math.cos(h * rad))

    # ---------- FIRST & LAST CONTACT (L1) ----------
    S: float = (o["a"] * o["v"] - o["u"] * o["b"]) / (o["n"] * o["L1p"])
    SS: float = math.sqrt(max(0.0, 1 - S * S))
    tau1: float = o["L1p"] / o["n"] * SS

    first_contact: float = t - tau1
    last_contact: float = t + tau1

    for _ in range(10):
        f = get_fundamental_arguments_unsafe(
            Xa,
            Ya,
            Da,
            Ma,
            L1a,
            L2a,
            tanf1,
            tanf2,
            first_contact,
            lat,
            lon,
            height,
            delta_t,
        )
        S: float = (f["a"] * f["v"] - f["u"] * f["b"]) / (f["n"] * f["L1p"])
        SS: float = math.sqrt(max(0.0, 1 - S * S))
        tf: float = (
            -(f["u"] * f["a"] + f["v"] * f["b"]) / (f["n"] ** 2)
            - f["L1p"] / f["n"] * SS
        )
        first_contact += tf

    for _ in range(10):
        f = get_fundamental_arguments_unsafe(
            Xa,
            Ya,
            Da,
            Ma,
            L1a,
            L2a,
            tanf1,
            tanf2,
            last_contact,
            lat,
            lon,
            height,
            delta_t,
        )
        S: float = (f["a"] * f["v"] - f["u"] * f["b"]) / (f["n"] * f["L1p"])
        SS: float = math.sqrt(max(0.0, 1 - S * S))
        tf: float = (
            -(f["u"] * f["a"] + f["v"] * f["b"]) / (f["n"] ** 2)
            + f["L1p"] / f["n"] * SS
        )
        last_contact += tf

    # ---------- SECOND & THIRD CONTACT (L2) ----------
    S: float = (o["a"] * o["v"] - o["u"] * o["b"]) / (o["n"] * o["L2p"])
    SS: float = math.sqrt(max(0.0, 1 - S * S))
    tau2: float = o["L2p"] / o["n"] * SS

    third_contact: float = t - tau2
    second_contact: float = t + tau2

    for _ in range(10):
        f = get_fundamental_arguments_unsafe(
            Xa,
            Ya,
            Da,
            Ma,
            L1a,
            L2a,
            tanf1,
            tanf2,
            third_contact,
            lat,
            lon,
            height,
            delta_t,
        )
        S: float = (f["a"] * f["v"] - f["u"] * f["b"]) / (f["n"] * f["L2p"])
        SS: float = math.sqrt(max(0.0, 1 - S * S))
        tf: float = (
            -(f["u"] * f["a"] + f["v"] * f["b"]) / (f["n"] ** 2)
            - f["L2p"] / f["n"] * SS
        )
        third_contact += tf

    for _ in range(10):
        f = get_fundamental_arguments_unsafe(
            Xa,
            Ya,
            Da,
            Ma,
            L1a,
            L2a,
            tanf1,
            tanf2,
            second_contact,
            lat,
            lon,
            height,
            delta_t,
        )
        S: float = (f["a"] * f["v"] - f["u"] * f["b"]) / (f["n"] * f["L2p"])
        SS: float = math.sqrt(max(0.0, 1 - S * S))
        tf: float = (
            -(f["u"] * f["a"] + f["v"] * f["b"]) / (f["n"] ** 2)
            + f["L2p"] / f["n"] * SS
        )
        second_contact += tf

    return {
        "t": t,
        "TTMaximum": T0 + t,
        "UT1Maximum": T0 + t - delta_t / 3600.0,
        "firstContact": first_contact,
        "secondContact": second_contact,
        "thirdContact": third_contact,
        "lastContact": last_contact,
        "magnitude": G,
        "ratio": A,
        "altitudeSunUnwrapped": h,
        "azimuthSunUnwrapped": q,
        "elements": o,
        "m": m,
    }


def get_local_circumstances(
    Xa, Ya, Da, Ma, L1a, L2a, tanf1, tanf2, T0, lat, lon, height, delta_t
):
    """
    Public wrapper that applies longitude sign convention
    and returns user-facing eclipse contact times.
    """
    e: dict = get_local_circumstances_unsafe(
        Xa, Ya, Da, Ma, L1a, L2a, tanf1, tanf2, T0, lat, -lon, height, delta_t
    )

    c1: float = None  # ty:ignore[invalid-assignment]
    c2: float = None  # ty:ignore[invalid-assignment]
    ge: float = None  # ty:ignore[invalid-assignment]
    c3: float = None  # ty:ignore[invalid-assignment]
    c4: float = None  # ty:ignore[invalid-assignment]
    mag: float = None  # ty:ignore[invalid-assignment]

    if e["magnitude"] > 0:
        c1: float = e["firstContact"]
        ge: float = e["t"]
        c4: float = e["lastContact"]
        mag: float = e["magnitude"]

        o: dict = e["elements"]
        if e["m"] < o["L2p"] or e["m"] < -o["L2p"]:
            c2: float = e["secondContact"]
            c3: float = e["thirdContact"]

    return c1, c2, ge, c3, c4, mag, e
