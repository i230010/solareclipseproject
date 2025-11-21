"""
main.py
-----------------

Demonstrates solar eclipse calculations, including eclipse search,
Besselian elements computation, maximum eclipse time, and gamma value.

"""

from skyfield.api import load
import pedatetime
import psefinder
import psebessel
import psegam


def main():
    """Main function to demonstrate eclipse computations and Besselian elements."""

    # Define the start and end datetime for eclipse search
    dt_start = pedatetime.datetime(2024, 1, 1, 0, 0, 0)
    dt_end = pedatetime.datetime(2025, 1, 1, 0, 0, 0)

    # =========================
    # Eclipse Search Example
    # =========================
    print("Eclipse Search Test:")
    psefinder.sefinder(dt_start, dt_end, False)
    # Uncomment below to print the separation angle as well
    # psefinder.sefinder(dt_start, dt_end, True)
    print("")

    # =========================
    # Besselian Elements Example
    # =========================
    dt_max = pedatetime.datetime(2024, 4, 8, 18, 17, 20)
    print(f"Data for {dt_max.isoformat()}Z Solar Eclipse")

    # Compute Besselian elements at surrounding times (2h before, 1h before, at max, etc.)
    TM2 = psebessel.besselian_find(dt_max - pedatetime.timedelta(0, 2, 0, 0))
    TM1 = psebessel.besselian_find(dt_max - pedatetime.timedelta(0, 1, 0, 0))
    T = psebessel.besselian_find(dt_max)
    TP1 = psebessel.besselian_find(dt_max + pedatetime.timedelta(0, 1, 0, 0))
    TP2 = psebessel.besselian_find(dt_max + pedatetime.timedelta(0, 2, 0, 0))

    # Generate Besselian polynomial coefficients
    Xarr = psebessel.find_besselian_polynomial(TM2[0], TM1[0], T[0], TP1[0], TP2[0])
    Yarr = psebessel.find_besselian_polynomial(TM2[1], TM1[1], T[1], TP1[1], TP2[1])
    Darr = psebessel.find_besselian_polynomial(TM2[2], TM1[2], T[2], TP1[2], TP2[2])
    L1arr = psebessel.find_besselian_polynomial(TM2[3], TM1[3], T[3], TP1[3], TP2[3])
    L2arr = psebessel.find_besselian_polynomial(TM2[4], TM1[4], T[4], TP1[4], TP2[4])
    Micro = psebessel.find_besselian_polynomial(TM2[5], TM1[5], T[5], TP1[5], TP2[5])

    # Unpack polynomial coefficients
    X0, X1, X2, X3 = Xarr
    Y0, Y1, Y2, Y3 = Yarr
    D0, D1, D2, D3 = Darr
    L10, L11, L12, L13 = L1arr
    L20, L21, L22, L23 = L2arr
    Micro0, Micro1, Micro2, Micro3 = Micro
    Tanf1, Tanf2 = T[6], T[7]

    # Print polynomial coefficients
    print(f"{'n':<3} {'X':>12} {'Y':>14} {'D':>14} {'L1':>14} {'L2':>14} {'Micro':>14}")
    print(
        f"0 {X0:14.10f} {Y0:14.10f} {D0:14.10f} {L10:14.10f} {L20:14.10f} {Micro0:14.10f}"
    )
    print(
        f"1 {X1:14.10f} {Y1:14.10f} {D1:14.10f} {L11:14.10f} {L21:14.10f} {Micro1:14.10f}"
    )
    print(
        f"2 {X2:14.10f} {Y2:14.10f} {D2:14.10f} {L12:14.10f} {L22:14.10f} {Micro2:14.10f}"
    )
    print(
        f"3 {X3:14.10f} {Y3:14.10f} {D3:14.10f} {L13:14.10f} {L23:14.10f} {Micro3:14.10f}"
    )
    print(f"tan(f1) = {Tanf1:14.10f} tan(f2) = {Tanf2:14.10f}")
    print("")

    # =========================
    # Maximum Eclipse Time
    # =========================
    ts = load.timescale()
    t = ts.utc(
        dt_max.year,
        dt_max.month,
        dt_max.day,
        dt_max.hour,
        dt_max.minute,
        dt_max.second,
    )

    t_hour = t.tt_strftime("%H")
    print(f"Maximum Eclipse: {t.tt_strftime()}")
    print(f"T0: {t_hour}h")
    print(f"Delta T: {t.delta_t}s")
    print("")

    # =========================
    # Gamma Value for Eclipse
    # =========================
    print(f"Gamma: {psegam.gamma(dt_max)}")
    print("")


if __name__ == "__main__":
    main()
