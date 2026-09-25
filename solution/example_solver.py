"""Example solver: probe the black box, then fit a hypothesised model.

Strategy (what an expert would do after exploratory probing):
  1. The output lives in ~5..11, moves by exactly +1 per decade of value_b/value_a
     in the mid range and is odd-symmetric in log(value_b/value_a) -> a buffer pH,
     pH ~ pK + log10(value_b / value_a).  value_a, value_b are acid/base concentrations.
  2. value_d (0..50) shifts the pK smoothly and strongly (~ -0.02 per unit): temperature
     in degC; pK(T) of an amine-type acid (large enthalpy of ionisation).
  3. value_c (inert salt) RAISES pH, roughly as sqrt(I)/(1+sqrt(I)) with a turnover
     at high I -> an activity correction on pH = -log a_H with a charge-symmetric
     acid BH+ <-> B + H+ (only gamma_H survives): Davies equation, A(T) growing with T.
  4. At extreme base excess the +1/decade law bends down -> water autoionisation
     in the charge balance.
  Model: Davies (0.3) with a Debye-Hueckel slope from eps_r(T) of water
  (Malmberg-Maryott), exact charge balance, and two Harned-Robinson laws
  pK(T) = a/T + b + cT fitted from probes (one for the buffer acid, one for water).

This solver is illustrative, not guaranteed optimal. Standard library only.
"""
import math

PROBE_BUDGET = 200

T_ZERO = 273.15
TEMPS = (0.0, 10.0, 20.0, 30.0, 40.0, 50.0)

_params = None


def _a_dh(t_c):
    # Debye-Hueckel slope (log10, molar basis) from water's permittivity (Malmberg & Maryott 1956)
    e, kb, na, eps0 = 1.602176634e-19, 1.380649e-23, 6.02214076e23, 8.8541878188e-12
    eps_r = 87.740 - 0.40008 * t_c + 9.398e-4 * t_c ** 2 - 1.410e-6 * t_c ** 3
    ekt = eps0 * eps_r * kb * (t_c + T_ZERO)
    return e ** 2 / (8 * math.pi * ekt) * math.sqrt(2 * na * 1000.0 * e ** 2 / ekt) / math.log(10)


def _model(c_bh, c_b, c_na, t_c, pk, pkw):
    a_dh, k, kw = _a_dh(t_c), 10 ** -pk, 10 ** -pkw
    c_tot, c_cl = c_bh + c_b, c_bh + c_na
    ionic, h, g = c_cl, 1e-8, 1.0
    for _ in range(40):
        s = math.sqrt(ionic)
        g = 10 ** (-a_dh * (s / (1 + s) - 0.3 * ionic))
        lo, hi = -14.0, -1.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            hh = 10 ** mid
            if hh + c_tot * hh / (hh + k) + c_na - c_cl - kw / (g * g * hh) > 0:
                hi = mid
            else:
                lo = mid
        h = 10 ** (0.5 * (lo + hi))
        ionic = c_cl + kw / (g * g * h)
    return -math.log10(g * h)


def _invert(target, f, lo, hi):
    """1-D bisection for x with f(x) = target, f monotone increasing."""
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _fit_hr(temps, pks):
    """Least-squares fit of pK = a/T + b + c*T (Harned-Robinson form), 3x3 normal equations."""
    rows = [(1.0 / (t + T_ZERO), 1.0, t + T_ZERO) for t in temps]
    m = [[sum(r[i] * r[j] for r in rows) for j in range(3)] for i in range(3)]
    v = [sum(r[i] * y for r, y in zip(rows, pks)) for i in range(3)]
    for i in range(3):  # Gaussian elimination
        for j in range(i + 1, 3):
            f = m[j][i] / m[i][i]
            m[j] = [mj - f * mi for mj, mi in zip(m[j], m[i])]
            v[j] -= f * v[i]
    x = [0.0, 0.0, 0.0]
    for i in (2, 1, 0):
        x[i] = (v[i] - sum(m[i][j] * x[j] for j in range(i + 1, 3))) / m[i][i]
    return x


def fit(query):
    global _params

    def q(a, b, c, d):
        return query({"value_a": a, "value_b": b, "value_c": c, "value_d": d})["output"]

    # Buffer pK(T): equimolar 0.05 M, no salt (water term negligible); a larger ratio
    # would couple in pK_w, so use the mid-point where d(pH)/d(pK) = 1.
    pk_w_guess = 14.0
    y_eq = [q(0.05, 0.05, 0.0, t) for t in TEMPS]
    pks = []
    for t, y in zip(TEMPS, y_eq):
        pks.append(_invert(y, lambda p, t=t: _model(0.05, 0.05, 0.0, t, p, pk_w_guess), 5.0, 11.0))
    hr_bh = _fit_hr(TEMPS, pks)

    # Water pK_w(T): maximal base excess and dilution, where [OH-] eats into BH+.
    pkws = []
    for t in TEMPS:
        pk = hr_bh[0] / (t + T_ZERO) + hr_bh[1] + hr_bh[2] * (t + T_ZERO)
        y = q(0.0005, 0.1, 0.0, t)
        # pH falls as K_w grows, i.e. pH increases with pK_w
        pkws.append(_invert(y, lambda w, t=t, pk=pk: _model(0.0005, 0.1, 0.0, t, pk, w), 12.0, 16.0))
    hr_w = _fit_hr(TEMPS, pkws)

    # Refine pK_bh once with the fitted water law (tiny effect at equimolar 0.05 M).
    pks = []
    for t, y in zip(TEMPS, y_eq):
        pkw = hr_w[0] / (t + T_ZERO) + hr_w[1] + hr_w[2] * (t + T_ZERO)
        pks.append(_invert(y, lambda p, t=t, w=pkw: _model(0.05, 0.05, 0.0, t, p, w), 5.0, 11.0))
    hr_bh = _fit_hr(TEMPS, pks)
    _params = (hr_bh, hr_w)


def solve(inputs):
    if _params is None:
        raise RuntimeError("call fit(query) first")
    hr_bh, hr_w = _params
    a, b, c, d = (float(inputs[k]) for k in ("value_a", "value_b", "value_c", "value_d"))
    temp = d + T_ZERO
    pk = hr_bh[0] / temp + hr_bh[1] + hr_bh[2] * temp
    pkw = hr_w[0] / temp + hr_w[1] + hr_w[2] * temp
    return {"output": round(_model(a, b, c, d, pk, pkw), 4)}
