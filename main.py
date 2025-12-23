"""
main.py
-------

Demonstrates solar eclipse calculations, including:

- Eclipse search (optional / disabled)
- Besselian elements computation
- Maximum eclipse time and gamma value
- Geographic coordinates of maximum eclipse
- Eclipse start and end times
- Plotting the eclipse path
"""

from typing import List

from skyfield.api import GREGORIAN_START, load

import pedatetime
import psebessel
import psegam
import psecentralcoords
import psestartendtime


def decimal_hours(hours: int, minutes: int, seconds: int) -> float:
    """
    Convert hours, minutes, and seconds to decimal hours.

    Example:
        0h 1m 30s -> 0.025 hours
    """
    return hours + minutes / 60.0 + seconds / 3600.0


def round2(n: float) -> int:
    """
    Standard rounding

    Example:
        3.4 -> 3
    """
    if n >= 0:
        return int(n + 0.5)
    else:
        return int(n - 0.5)


def main() -> None:
    """
    Main driver function demonstrating eclipse computations
    and visualization.
    """

    # ------------------------------------------------------------------
    # Eclipse Search Example (disabled)
    # ------------------------------------------------------------------
    # Uncomment to perform an actual eclipse search.
    #
    # print("Solar Eclipse Search:")
    # import psefinder
    # dt_start = pedatetime.datetime(2024, 1, 1, 0, 0, 0)
    # dt_end = pedatetime.datetime(2025, 1, 1, 0, 0, 0)
    #
    # step = pedatetime.timedelta(0, 2, 0, 0)
    #
    # psefinder.sefinder(dt_start, dt_end, step, False, False)
    # psefinder.sefinder(dt_start, dt_end, step, False, True)  # include separation angle in radians
    # psefinder.sefinder(dt_start, dt_end, step, True, False)  # TT
    # psefinder.sefinder(dt_start, dt_end, step, True, True)  # sep angle + TT

    print()

    # ------------------------------------------------------------------
    # Besselian Elements Example
    # ------------------------------------------------------------------
    ts = load.timescale()
    ts.julian_calendar_cutoff = GREGORIAN_START

    # UT1 time of maximum eclipse
    dt_max_ut1 = pedatetime.datetime(2024, 4, 8, 18, 17, 20)
    t_max = ts.ut1(
        dt_max_ut1.year,
        dt_max_ut1.month,
        dt_max_ut1.day,
        dt_max_ut1.hour,
        dt_max_ut1.minute,
        dt_max_ut1.second,
    )

    # Delta-T in seconds (TT - UT1)
    delta_t_sec: int = round2(float(t_max.delta_t))
    dt_max_tt = dt_max_ut1 + pedatetime.timedelta(0, 0, 0, delta_t_sec)

    print(f"Data for {dt_max_ut1.isoformat()} UT1 ({dt_max_tt.isoformat()} TT) Eclipse")

    # ------------------------------------------------------------------
    # Round to nearest hour for polynomial interpolation
    # ------------------------------------------------------------------
    dt_max_rounded = dt_max_tt.copy()

    if dt_max_tt.minute >= 30:
        dt_max_rounded.sub_minutes(dt_max_tt.minute)
        dt_max_rounded.sub_seconds(dt_max_tt.second)
        dt_max_rounded.add_hour()
    else:
        dt_max_rounded.sub_minutes(dt_max_tt.minute)
        dt_max_rounded.sub_seconds(dt_max_tt.second)

    # ------------------------------------------------------------------
    # Compute Besselian elements at plus/minus 2 hours
    # ------------------------------------------------------------------
    TM2 = psebessel.besselian_find(dt_max_rounded - pedatetime.timedelta(0, 2, 0, 0))
    TM1 = psebessel.besselian_find(dt_max_rounded - pedatetime.timedelta(0, 1, 0, 0))
    T0 = psebessel.besselian_find(dt_max_rounded)
    TP1 = psebessel.besselian_find(dt_max_rounded + pedatetime.timedelta(0, 1, 0, 0))
    TP2 = psebessel.besselian_find(dt_max_rounded + pedatetime.timedelta(0, 2, 0, 0))

    # ------------------------------------------------------------------
    # Polynomial coefficients (Besselian elements)
    # For L1 and L2, use a second-degree polynomial. However for improved accuracy, use a third-degree polynomial
    # ------------------------------------------------------------------
    X_poly = psebessel.find_3rd_degree_polynomial(TM2[0], TM1[0], T0[0], TP1[0], TP2[0])
    Y_poly = psebessel.find_3rd_degree_polynomial(TM2[1], TM1[1], T0[1], TP1[1], TP2[1])
    D_poly = psebessel.find_3rd_degree_polynomial(TM2[2], TM1[2], T0[2], TP1[2], TP2[2])
    L1_poly = psebessel.find_3rd_degree_polynomial(
        TM2[3], TM1[3], T0[3], TP1[3], TP2[3]
    )
    L2_poly = psebessel.find_3rd_degree_polynomial(
        TM2[4], TM1[4], T0[4], TP1[4], TP2[4]
    )

    # Micro is linear
    Micro_poly = psebessel.find_1st_degree_polynomial(T0[5], TP1[5])

    tan_f1: float = T0[6]
    tan_f2: float = T0[7]

    # ------------------------------------------------------------------
    # Display coefficients
    # ------------------------------------------------------------------
    print(f"{'n':<3} {'X':>12} {'Y':>14} {'D':>14} {'L1':>14} {'L2':>14} {'Micro':>14}")

    for n, values in enumerate(
        zip(X_poly, Y_poly, D_poly, L1_poly, L2_poly, Micro_poly)
    ):
        print(
            f"{n} "
            f"{values[0]:14.10f} {values[1]:14.10f} {values[2]:14.10f} "
            f"{values[3]:14.10f} {values[4]:14.10f} {values[5]:14.10f}"
        )

    print(f"tan(f1) = {tan_f1:14.10f}  tan(f2) = {tan_f2:14.10f}")

    # ------------------------------------------------------------------
    # Maximum Eclipse Time
    # ------------------------------------------------------------------
    decimal_time_tt: float = (dt_max_tt - dt_max_rounded.copy()).total_seconds / 3600.0

    print(f"Maximum Eclipse: {dt_max_tt.isoformat()} TT")
    print(f"Delta T: {delta_t_sec}s")

    # ------------------------------------------------------------------
    # Gamma at Maximum Eclipse
    # ------------------------------------------------------------------
    gamma_val: float = psegam.gamma(X_poly, Y_poly, decimal_time_tt)
    print(f"Gamma: {gamma_val}")

    # ------------------------------------------------------------------
    # Maximum Eclipse Geographic Location
    # ------------------------------------------------------------------
    lat_max, lon_max = psecentralcoords.coords(
        X_poly,
        Y_poly,
        D_poly,
        Micro_poly,
        delta_t_sec,
        decimal_time_tt,
    )

    if lat_max is not None and lon_max is not None:
        print(f"Maximum Eclipse Location: {lat_max}, {lon_max}")

    # ------------------------------------------------------------------
    # Eclipse Start and End Times
    # ------------------------------------------------------------------
    base_dt_hour = dt_max_rounded.copy()
    print(f"T0: {base_dt_hour.isoformat()} TT")

    # Penumbral contacts
    t_start_pen, t_end_pen = psestartendtime.startendtime(X_poly, Y_poly, L1_poly)

    tt_start_pen = base_dt_hour + pedatetime.timedelta(
        0, 0, 0, int(round2(t_start_pen * 3600))
    )
    tt_end_pen = base_dt_hour + pedatetime.timedelta(
        0, 0, 0, int(round2(t_end_pen * 3600))
    )

    print(f"Eclipse Start (Penumbra): {tt_start_pen.isoformat()} TT")
    print(f"Eclipse End   (Penumbra): {tt_end_pen.isoformat()} TT")

    # Umbral contacts (only if central eclipse exists)
    if lat_max is not None and lon_max is not None:
        t_start_umb, t_end_umb = psestartendtime.startendtime(X_poly, Y_poly, L2_poly)

        tt_start_umb = base_dt_hour + pedatetime.timedelta(
            0, 0, 0, int(round2(t_start_umb * 3600))
        )
        tt_end_umb = base_dt_hour + pedatetime.timedelta(
            0, 0, 0, int(round2(t_end_umb * 3600))
        )

        print(f"Eclipse Start (Central Umbra): {tt_start_umb.isoformat()} TT")
        print(f"Eclipse End   (Central Umbra): {tt_end_umb.isoformat()} TT")

    # ------------------------------------------------------------------
    # Compute Eclipse Path for Plotting
    # ------------------------------------------------------------------
    if lat_max is not None and lon_max is not None:
        step_seconds: int = 60
        step_hours: float = decimal_hours(0, 0, step_seconds)

        path_lats: List[float] = []
        path_lons: List[float] = []

        current_time: float = t_start_pen

        while current_time < t_end_pen:
            current_time += step_hours

            lat, lon = psecentralcoords.coords(
                X_poly,
                Y_poly,
                D_poly,
                Micro_poly,
                delta_t_sec,
                current_time,
            )

            if lat is not None and lon is not None:
                path_lats.append(lat)
                path_lons.append(lon)

        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature

        proj = ccrs.Orthographic(central_longitude=lon_max, central_latitude=lat_max)

        fig = plt.figure(figsize=(8, 8))
        ax = plt.axes(projection=proj)

        ax.add_feature(cfeature.OCEAN, facecolor="blue")
        ax.add_feature(cfeature.LAND, facecolor="green")
        ax.add_feature(cfeature.LAKES, facecolor="lightblue")
        ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=1)
        ax.coastlines(color="black", linewidth=0.8)

        ax.plot(
            path_lons,
            path_lats,
            color="black",
            linewidth=2,
            transform=ccrs.PlateCarree(),
        )

        ax.plot(
            lon_max,
            lat_max,
            marker="o",
            color="white",
            markersize=6,
            transform=ccrs.PlateCarree(),
        )

        ax.set_global()
        fig.savefig("central_path.png")


if __name__ == "__main__":
    main()
