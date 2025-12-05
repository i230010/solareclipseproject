"""
main.py
-------

Demonstrates solar eclipse calculations, including eclipse search,
Besselian elements computation, maximum eclipse time, gamma value, and
plotting the eclipse path.
"""

from skyfield.api import load

import pedatetime
import psebessel
import psegam
import pcentercoords
import pstartendtimel
import psecentral


def decimal_hours(hours: int, minutes: int, seconds: int) -> float:
    """Convert hours, minutes, and seconds to decimal hours."""
    return hours + minutes / 60 + seconds / 3600


def main():
    """Main function demonstrating eclipse computations and visualization."""

    # ------------------------------------------------------------------
    # Eclipse Search Example (disabled)
    # ------------------------------------------------------------------
    # Uncomment to perform an actual eclipse search.
    #
    # print("Eclipse Search Example:")
    # import psefinder
    #
    # dt_start = pedatetime.datetime(2020, 1, 1, 0, 0, 0)
    # dt_end = pedatetime.datetime(2020, 7, 1, 0, 0, 0)
    #
    # psefinder.sefinder(dt_start, dt_end, False, pedatetime.timedelta(0, 2, 0, 0))
    # psefinder.sefinder(dt_start, dt_end, True, pedatetime.timedelta(0, 2, 0, 0))  # includes sep angle

    print("")

    # ------------------------------------------------------------------
    # Besselian Elements Example
    # ------------------------------------------------------------------
    # dt_max = pedatetime.datetime(2024, 4, 8, 18, 17, 20)
    dt_max = pedatetime.datetime(2020, 6, 21, 6, 40, 6)
    print(f"Data for {dt_max.isoformat()}Z Solar Eclipse")

    # Round down to the nearest hour for polynomial interpolation
    dt_max_rounded = dt_max.copy()
    if dt_max.minute >= 30:
        dt_max_rounded.sub_minutes(dt_max.minute)
        dt_max_rounded.sub_seconds(dt_max.second)
        dt_max_rounded.add_hour()
    else :
        dt_max_rounded.sub_minutes(dt_max.minute)
        dt_max_rounded.sub_seconds(dt_max.second)
    
    # Compute Besselian elements at times -2h, -1h, 0h, +1h, +2h
    TM2 = psebessel.besselian_find(
        dt_max_rounded - pedatetime.timedelta(0, 2, 0, 0)
    )
    TM1 = psebessel.besselian_find(
        dt_max_rounded - pedatetime.timedelta(0, 1, 0, 0)
    )
    T = psebessel.besselian_find(dt_max_rounded)
    TP1 = psebessel.besselian_find(
        dt_max_rounded + pedatetime.timedelta(0, 1, 0, 0)
    )
    TP2 = psebessel.besselian_find(
        dt_max_rounded + pedatetime.timedelta(0, 2, 0, 0)
    )

    # Generate polynomial coefficients
    X_poly = psebessel.find_besselian_polynomial(
        TM2[0], TM1[0], T[0], TP1[0], TP2[0]
    )
    Y_poly = psebessel.find_besselian_polynomial(
        TM2[1], TM1[1], T[1], TP1[1], TP2[1]
    )
    D_poly = psebessel.find_besselian_polynomial(
        TM2[2], TM1[2], T[2], TP1[2], TP2[2]
    )
    L1_poly = psebessel.find_besselian_polynomial(
        TM2[3], TM1[3], T[3], TP1[3], TP2[3]
    )
    L2_poly = psebessel.find_besselian_polynomial(
        TM2[4], TM1[4], T[4], TP1[4], TP2[4]
    )
    Micro_poly = psebessel.find_besselian_polynomial(
        TM2[5], TM1[5], T[5], TP1[5], TP2[5]
    )

    # Unpack polynomial coefficients
    tan_f1, tan_f2 = T[6], T[7]

    # Display coefficients
    print(
        f"{'n':<3} {'X':>12} {'Y':>14} {'D':>14} "
        f"{'L1':>14} {'L2':>14} {'Micro':>14}"
    )
    for n, vals in enumerate(
        zip(X_poly, Y_poly, D_poly, L1_poly, L2_poly, Micro_poly)
    ):
        print(
            f"{n} "
            f"{vals[0]:14.10f} {vals[1]:14.10f} {vals[2]:14.10f} "
            f"{vals[3]:14.10f} {vals[4]:14.10f} {vals[5]:14.10f}"
        )

    print(f"tan(f1) = {tan_f1:14.10f}  tan(f2) = {tan_f2:14.10f}")

    # ------------------------------------------------------------------
    # Maximum Eclipse Time
    # ------------------------------------------------------------------
    ts = load.timescale()
    t_max = ts.utc(
        dt_max.year, dt_max.month, dt_max.day,
        dt_max.hour, dt_max.minute, dt_max.second
    )
    
    #T0
    decimal_time = (
        decimal_hours(
            int(t_max.tt_strftime("%H")),
            int(t_max.tt_strftime("%M")),
            int(t_max.tt_strftime("%S")),
        )
        - dt_max_rounded.hour
    )

    delta_t = t_max.delta_t
    print(f"Maximum Eclipse: {t_max.tt_strftime()}")
    print(f"Delta T: {delta_t}s")

    # ------------------------------------------------------------------
    # Gamma at Maximum Eclipse
    # ------------------------------------------------------------------

    gamma_val = psegam.gamma(X_poly, Y_poly, decimal_time)
    print(f"Gamma: {gamma_val}")

    print(
        "Central: "
        f"{psecentral.central(gamma_val)}, "
        "Umbra Touches Earth: "
        f"{psecentral.umbra_touch(gamma_val)}, "
        "Eclipse Exists: "
        f"{psecentral.exist(gamma_val)}"
    )

    # ------------------------------------------------------------------
    # Maximum Eclipse Geographic Location
    # ------------------------------------------------------------------
    lat_max, lon_max = pcentercoords.coords(
        X_poly, Y_poly, D_poly, Micro_poly, delta_t, decimal_time
    )
    print(f"Maximum Eclipse Location: {lat_max}, {lon_max}")

    # ------------------------------------------------------------------
    # Eclipse Start and End Times
    # ------------------------------------------------------------------
    base_dt_hour = dt_max_rounded.copy()

    print(f"T0: {base_dt_hour.isoformat()} TT")

    # Penumbra
    t_start_pen, t_end_pen = pstartendtimel.startendtime(
        X_poly, Y_poly, L1_poly
    )

    tt_start_pen = base_dt_hour + pedatetime.timedelta(
        0, 0, 0, int(t_start_pen * 3600)
    )
    tt_end_pen = base_dt_hour + pedatetime.timedelta(
        0, 0, 0, int(t_end_pen * 3600)
    )

    print(f"Eclipse Start (Penumbra): {tt_start_pen.isoformat()} TT")
    print(f"Eclipse End (Penumbra):   {tt_end_pen.isoformat()} TT")

    # Umbra
    t_start_umb, t_end_umb = pstartendtimel.startendtime(
        X_poly, Y_poly, L2_poly
    )

    tt_start_umb = base_dt_hour + pedatetime.timedelta(
        0, 0, 0, int(t_start_umb * 3600)
    )
    tt_end_umb = base_dt_hour + pedatetime.timedelta(
        0, 0, 0, int(t_end_umb * 3600)
    )

    print(f"Eclipse Start (Umbra):    {tt_start_umb.isoformat()} TT")
    print(f"Eclipse End (Umbra):      {tt_end_umb.isoformat()} TT")

    # ------------------------------------------------------------------
    # Compute Eclipse Path for Plotting
    # ------------------------------------------------------------------
    import plotly.graph_objects as go
    step_sec = 60
    start_time = t_start_pen
    end_time = t_end_pen

    path_lats = []
    path_lons = []

    current_time = start_time
    step_hours = decimal_hours(0, 0, step_sec)

    while current_time < end_time:
        current_time += step_hours
        lat, lon = pcentercoords.coords(
            X_poly, Y_poly, D_poly, Micro_poly, delta_t, current_time
        )
        if lat is not None and lon is not None:
            path_lats.append(lat)
            path_lons.append(lon)
            
#     ------------------------------------------------------------------
#     Plot Eclipse Path
#     ------------------------------------------------------------------
    fig = go.Figure()

    fig.add_trace(
        go.Scattergeo(
            lat=path_lats,
            lon=path_lons,
            mode="lines",
            line=dict(width=2, color="black"),
        )
    )

    fig.update_geos(
        projection_type="orthographic",
        projection_rotation=dict(lat=lat_max, lon=lon_max),
        projection_scale=0.80,
        showland=True,
        landcolor="green",
        showocean=True,
        oceancolor="skyblue",
        showcountries=True,
        countrycolor="black",
        showlakes=True,
        lakecolor="blue",
    )

    fig.update_layout(
        width=1000,
        height=1000,
        margin=dict(l=0, r=0, t=0, b=0),
        title="Eclipse Path",
    )

    fig.write_image("path.png")


if __name__ == "__main__":
    main()
