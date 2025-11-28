"""
main.py
-----------------

Demonstrates solar eclipse calculations, including eclipse search,
Besselian elements computation, maximum eclipse time, gamma value,
and plotting the eclipse path.
"""

from skyfield.api import load
import pedatetime
import psebessel
import psegam
import pcentercoords
import pstartendtimel
import plotly.graph_objects as go


def decimal_hours(hours: int, minutes: int, seconds: int) -> float:
    """Convert hours, minutes, and seconds to decimal hours."""
    return hours + minutes / 60 + seconds / 3600


def main():
    """Main function demonstrating eclipse computations and visualization."""

    # Uncomment to perform actual eclipse search
    """
    print("Eclipse Search Example:")
    import psefinder

    # ------------------------
    # Define Eclipse Search Period
    # ------------------------

    dt_start = pedatetime.datetime(2024, 1, 1, 0, 0, 0)
    dt_end = pedatetime.datetime(2025, 1, 1, 0, 0, 0)
    
    psefinder.sefinder(dt_start, dt_end, False, 2)
    # psefinder.sefinder(dt_start, dt_end, True, 2)  # includes separation angle
    """

    print("")

    # ------------------------
    # Besselian Elements Example
    # ------------------------
    dt_max = pedatetime.datetime(2024, 4, 8, 18, 17, 20)
    print(f"Data for {dt_max.isoformat()}Z Solar Eclipse")

    # Round down to the nearest hour for polynomial interpolation
    dt_max_rounded = dt_max.copy()
    dt_max_rounded.sub_minutes(dt_max.minute)
    dt_max_rounded.sub_seconds(dt_max.second)

    # Compute Besselian elements at surrounding times (-2h, -1h, max, +1h, +2h)
    TM2 = psebessel.besselian_find(dt_max_rounded - pedatetime.timedelta(0, 2, 0, 0))
    TM1 = psebessel.besselian_find(dt_max_rounded - pedatetime.timedelta(0, 1, 0, 0))
    T = psebessel.besselian_find(dt_max_rounded)
    TP1 = psebessel.besselian_find(dt_max_rounded + pedatetime.timedelta(0, 1, 0, 0))
    TP2 = psebessel.besselian_find(dt_max_rounded + pedatetime.timedelta(0, 2, 0, 0))

    # Generate Besselian polynomial coefficients
    X_poly = psebessel.find_besselian_polynomial(TM2[0], TM1[0], T[0], TP1[0], TP2[0])
    Y_poly = psebessel.find_besselian_polynomial(TM2[1], TM1[1], T[1], TP1[1], TP2[1])
    D_poly = psebessel.find_besselian_polynomial(TM2[2], TM1[2], T[2], TP1[2], TP2[2])
    L1_poly = psebessel.find_besselian_polynomial(TM2[3], TM1[3], T[3], TP1[3], TP2[3])
    L2_poly = psebessel.find_besselian_polynomial(TM2[4], TM1[4], T[4], TP1[4], TP2[4])
    Micro_poly = psebessel.find_besselian_polynomial(
        TM2[5], TM1[5], T[5], TP1[5], TP2[5]
    )

    # Unpack polynomial coefficients for display
    X0, X1, X2, X3 = X_poly
    Y0, Y1, Y2, Y3 = Y_poly
    D0, D1, D2, D3 = D_poly
    L10, L11, L12, L13 = L1_poly
    L20, L21, L22, L23 = L2_poly
    Micro0, Micro1, Micro2, Micro3 = Micro_poly
    tan_f1, tan_f2 = T[6], T[7]

    # Print polynomial coefficients
    print(f"{'n':<3} {'X':>12} {'Y':>14} {'D':>14} {'L1':>14} {'L2':>14} {'Micro':>14}")
    for n, vals in enumerate(zip(X_poly, Y_poly, D_poly, L1_poly, L2_poly, Micro_poly)):
        print(
            f"{n} {vals[0]:14.10f} {vals[1]:14.10f} {vals[2]:14.10f} {vals[3]:14.10f} {vals[4]:14.10f} {vals[5]:14.10f}"
        )
    print(f"tan(f1) = {tan_f1:14.10f} tan(f2) = {tan_f2:14.10f}")

    # ------------------------
    # Maximum Eclipse Time
    # ------------------------
    ts = load.timescale()
    t_max = ts.utc(
        dt_max.year, dt_max.month, dt_max.day, dt_max.hour, dt_max.minute, dt_max.second
    )

    delta_t = t_max.delta_t  # ΔT in seconds
    print(f"Maximum Eclipse: {t_max.tt_strftime()}")
    print(f"Delta T: {delta_t}s")

    # ------------------------
    # Gamma Value at Maximum Eclipse
    # ------------------------
    decimal_time = (
        decimal_hours(
            int(t_max.tt_strftime("%H")),
            int(t_max.tt_strftime("%M")),
            int(t_max.tt_strftime("%S")),
        )
        - dt_max.hour
    )
    gamma_val = psegam.gamma(X_poly, Y_poly, decimal_time)
    print(f"Gamma: {gamma_val}")

    # Compute geographic location of maximum eclipse
    lat_max, lon_max = pcentercoords.coords(
        X_poly, Y_poly, D_poly, Micro_poly, delta_t, decimal_time
    )
    print(f"Maximum Eclipse Location: {lat_max}, {lon_max}")

    # ------------------------
    # Eclipse Start and End Times
    # ------------------------
    base_dt_hour = pedatetime.datetime(
        dt_max.year, dt_max.month, dt_max.day, dt_max.hour, 0, 0
    )

    # Penumbra
    t_start_penumbra, t_end_penumbra = pstartendtimel.startendtime(
        X_poly, Y_poly, L1_poly
    )
    tt_start_penumbra = base_dt_hour + pedatetime.timedelta(
        0, 0, 0, int(t_start_penumbra * 3600)
    )
    tt_end_penumbra = base_dt_hour + pedatetime.timedelta(
        0, 0, 0, int(t_end_penumbra * 3600)
    )
    print(f"Eclipse Start (Penumbra): {tt_start_penumbra.isoformat()} TT")
    print(f"Eclipse End (Penumbra): {tt_end_penumbra.isoformat()} TT")

    # Umbra
    t_start_umbra, t_end_umbra = pstartendtimel.startendtime(X_poly, Y_poly, L2_poly)
    tt_start_umbra = base_dt_hour + pedatetime.timedelta(
        0, 0, 0, int(t_start_umbra * 3600)
    )
    tt_end_umbra = base_dt_hour + pedatetime.timedelta(0, 0, 0, int(t_end_umbra * 3600))
    print(f"Eclipse Start (Umbra): {tt_start_umbra.isoformat()} TT")
    print(f"Eclipse End (Umbra): {tt_end_umbra.isoformat()} TT")

    # ------------------------
    # Compute Eclipse Path Coordinates
    # ------------------------
    step_sec = 1  # 1-second time step in decimal hours
    start_time = t_start_penumbra
    end_time = t_end_penumbra
    path_lats = []
    path_lons = []

    current_time = start_time
    while current_time < end_time:
        current_time += decimal_hours(0, 0, step_sec)
        lat, lon = pcentercoords.coords(
            X_poly, Y_poly, D_poly, Micro_poly, delta_t, current_time
        )
        if lat is not None and lon is not None:
            path_lats.append(lat)
            path_lons.append(lon)

    # ------------------------
    # Plot Eclipse Path
    # ------------------------
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
        width=3000, height=3000, margin=dict(l=0, r=0, t=0, b=0), title="Eclipse Path"
    )

    fig.write_image("eclipse_path.png", scale=3)


if __name__ == "__main__":
    main()
