import math

def eval_cubic(a, t):
    return a[0] + a[1]*t + a[2]*t*t + a[3]*t*t*t

def eval_quad(a, t):
    return a[0] + a[1]*t + a[2]*t*t

def eval_linear(a, t):
    return a[0] + a[1]*t

def coords(Xa, Ya, Da, Ma, delta_t, T):

    # Step 1: time

    # Step 2: Besselian evaluation
    X = eval_cubic(Xa, T)
    Y = eval_cubic(Ya, T)
    d = math.radians(eval_quad(Da, T))
    m = math.radians(eval_linear(Ma, T))   # Greenwich hour angle in radians

    # Earth ellipsoid constants
    e2 = 0.006694385
    one_minus_f = 0.99664719

    # Step 3: ellipsoid corrections
    omega = 1.0 / math.sqrt(1 - e2 * (math.cos(d)**2))
    y1 = omega * Y
    b1 = omega * math.sin(d)
    b2 = one_minus_f * omega * math.cos(d)

    # Step 4: compute B
    Bsq = (1 - X*X - y1*y1)
    if Bsq < 0:
        return None, None
    B = math.sqrt(Bsq)

    # Step 5: latitude components
    sinphi1 = B*b1 + y1*b2
    phi1 = math.asin(sinphi1)  # geocentric latitude

    # Oblateness correction
    phi = math.atan(1.00336409 * math.tan(phi1))  # geodetic latitude

    # Step 6: hour angle
    sinH = X / math.cos(phi1)
    cosH = (B*b2 - y1*b1) / math.cos(phi1)
    H = math.atan2(sinH, cosH)

    # Step 7: longitude
    lambda_geo = m - H - (0.00417807 * delta_t * math.pi/180)

    # Convert to degrees
    return math.degrees(phi),((-math.degrees(lambda_geo) + 180) % 360) - 180