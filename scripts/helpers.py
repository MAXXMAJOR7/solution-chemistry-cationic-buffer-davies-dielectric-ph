"""Shared utilities for task scripts and grader."""
import importlib.util
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GOLDEN_PATH = os.path.join(ROOT, "golden", "test_data.json")
ORACLE_PATH = os.path.join(ROOT, "oracle", "implement.py")

INPUT_KEYS = ("value_a", "value_b", "value_c", "value_d")
INPUT_RANGES = {
    "value_a": (0.0005, 0.1),
    "value_b": (0.0005, 0.1),
    "value_c": (0.0, 0.2),
    "value_d": (0.0, 50.0),
}
TOLERANCE = 0.001  # absolute, pH units (output spans ~5.1 .. 10.6)


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_oracle():
    return load_module(ORACLE_PATH, "oracle_impl")


def load_golden(path=GOLDEN_PATH):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def variant(c_bh, c_b, c_na, t_c, activity=True, charge_aware=True, thermal=True, full_balance=True):
    """Oracle with individual twists frozen, for twist-integrity checks and baselines.

    activity=False     -> ideal solution, gamma = 1 (textbook Henderson-Hasselbalch family)
    charge_aware=False -> treats the acid as neutral HA -> A- (acetate-like): the activity
                          correction enters as K_a' = K_a / gamma^2, i.e. the salt LOWERS pH
    thermal=False      -> every temperature law frozen at 25 degC (pK_bh, pK_w and A)
    full_balance=False -> neglect h and [OH-] in the charge balance (h = K_bh [BH+]/[B])
    """
    o = load_oracle()
    t_eval = t_c if thermal else 25.0
    temp = t_eval + o.T_ZERO
    k_bh = 10 ** -(o.BH_A / temp + o.BH_B - o.BH_C * temp)
    k_w = 10 ** (-o.W_A / temp + o.W_B - o.W_C * temp)
    a_dh = o.debye_hueckel_a(t_eval)
    c_tot, c_cl = c_bh + c_b, c_bh + c_na
    ionic, h, gamma = c_cl, None, 1.0
    for _ in range(o.N_FIXED_POINT):
        gamma = o.davies_gamma(a_dh, ionic) if activity else 1.0
        # effective acid constant on the concentration scale (h = [H+])
        k_eff = k_bh if charge_aware else k_bh / gamma ** 2
        if not full_balance:
            h = k_eff * c_bh / c_b
        else:
            def residual(log_h, g=gamma, k=k_eff):
                hh = 10 ** log_h
                return hh + c_tot * hh / (hh + k) + c_na - c_cl - k_w / (g ** 2 * hh)
            lo, hi = -14.0, -1.0
            for _ in range(o.N_BISECT):
                mid = 0.5 * (lo + hi)
                if residual(mid) > 0:
                    hi = mid
                else:
                    lo = mid
            h = 10 ** (0.5 * (lo + hi))
        ionic = c_cl + (k_w / (gamma ** 2 * h) if full_balance else 0.0)
    return round(-math.log10(gamma * h), 4)
