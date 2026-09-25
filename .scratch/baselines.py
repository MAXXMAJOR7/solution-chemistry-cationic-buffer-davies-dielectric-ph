"""Baseline models scored on golden data (not shipped)."""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from helpers import TOLERANCE, load_golden, variant  # noqa: E402

cases = load_golden()["cases"]


def hh25(a, b, c, d):
    return round(7.762 + math.log10(b / a), 4)


def vanthoff(a, b, c, d):
    # 25 degC pK with a constant-enthalpy van't Hoff shift (dH from dpK/dT at 25 degC), activities on
    t = d + 273.15
    dh_over_r = 1341.16 + 0.0045666 * 298.15 ** 2  # ln10 cancels in log form
    pk = 7.762 + dh_over_r * (1 / t - 1 / 298.15)
    return variant(a, b, c, 25.0) - 7.7619 + pk


models = {
    "Henderson-Hasselbalch, 25 degC pK": hh25,
    "Ideal (gamma = 1), exact T laws": lambda *x: variant(*x, activity=False),
    "Acetate-type activity (salt lowers pH)": lambda *x: variant(*x, charge_aware=False),
    "T laws frozen at 25 degC": lambda *x: variant(*x, thermal=False),
    "No water/h in charge balance": lambda *x: variant(*x, full_balance=False),
    "Constant-dH van't Hoff pK(T)": vanthoff,
}
for name, f in models.items():
    errs = sorted(abs(f(*(c["input"][k] for k in ("value_a", "value_b", "value_c", "value_d")))
                      - c["expected"]["output"]) for c in cases)
    ok = sum(e <= TOLERANCE for e in errs)
    print(f"{name:42s} {ok:2d}/{len(cases)}  median {errs[len(errs)//2]:.4f}  max {errs[-1]:.4f}")
