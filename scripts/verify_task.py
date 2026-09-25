"""Automated GATE 2 checks. Exits non-zero on any failure."""
import ast
import json
import math
import random
import sys

from helpers import (GOLDEN_PATH, INPUT_KEYS, INPUT_RANGES, ORACLE_PATH, load_golden,
                     load_oracle, variant)

# Every numeric literal allowed in oracle/implement.py, with its justification.
ALLOWED_CONSTANTS = {
    1.602176634e-19: "elementary charge (exact, SI 2019)",
    1.380649e-23: "Boltzmann constant (exact, SI 2019)",
    6.02214076e23: "Avogadro constant (exact, SI 2019)",
    8.8541878188e-12: "vacuum permittivity (CODATA 2022)",
    273.15: "0 degC in K (exact)",
    1000.0: "L per m^3 (exact)",
    1341.16: "Bates & Allen 1960 eq (5)", 4.6252: "Bates & Allen 1960 eq (5)",
    0.0045666: "Bates & Allen 1960 eq (5)",
    4470.99: "Harned & Robinson 1940 log K_w", 6.0875: "Harned & Robinson 1940 log K_w",
    0.01706: "Harned & Robinson 1940 log K_w",
    87.740: "Malmberg & Maryott 1956", 0.40008: "Malmberg & Maryott 1956",
    9.398e-4: "Malmberg & Maryott 1956", 1.410e-6: "Malmberg & Maryott 1956",
    0.3: "Davies 1962 linear-term coefficient",
    8: "Debye-Hueckel: 8 pi eps kT denominator (exact theory)",
    2: "structural: kappa^2 = 2 N_A e^2 I / eps kT; squares; Davies z^2 = 1",
    3: "structural: cubic term of eps_r(t)",
    1: "structural: 1 + sqrt(I)",
    10: "log10 base",
    0: "structural (domain check)",
    0.5: "bisection midpoint",
    50: "upper temperature of Bates & Allen fit (degC)",
    -14.0: "bisection bracket, log10 h (pH 14: above any pH in the domain)",
    -1.0: "bisection bracket, log10 h (pH 1: below any pH in the domain)",
    200: "fixed bisection iteration count (determinism)",
    60: "fixed ionic-strength fixed-point iteration count (determinism)",
    4: "4-decimal output rounding (task spec)",
}

failures = []


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not cond:
        failures.append(name)


def numeric_literals(path):
    tree = ast.parse(open(path, encoding="utf-8").read())
    vals = []
    for node in ast.walk(tree):
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub) and isinstance(node.operand, ast.Constant):
            if isinstance(node.operand.value, (int, float)) and not isinstance(node.operand.value, bool):
                vals.append(-node.operand.value)
                node.operand.value = None  # avoid double counting the positive literal
        elif isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
            vals.append(node.value)
    return vals


def sample(rng):
    lo, hi = (math.log10(x) for x in INPUT_RANGES["value_a"])
    return (10 ** rng.uniform(lo, hi), 10 ** rng.uniform(lo, hi),
            rng.uniform(*INPUT_RANGES["value_c"]), rng.uniform(*INPUT_RANGES["value_d"]))


def main():
    o = load_oracle()
    rng = random.Random(7)
    pts = [sample(rng) for _ in range(300)]

    # Determinism
    check("deterministic", all(o.compute(*p) == o.compute(*p) for p in pts))

    # Output schema
    out = o.oracle({"value_a": 0.01, "value_b": 0.01, "value_c": 0.1, "value_d": 25.0})
    check("output key is 'output' and scalar float",
          list(out) == ["output"] and isinstance(out["output"], float))
    check("4-decimal precision", all(round(o.compute(*p), 4) == o.compute(*p) for p in pts))

    # No invented constants
    unknown = [v for v in numeric_literals(ORACLE_PATH) if v not in ALLOWED_CONSTANTS]
    check("every oracle numeric literal is justified", not unknown, f"unjustified: {unknown}" if unknown else "")

    # Literature anchors
    t25 = 25 + o.T_ZERO
    pk25 = o.BH_A / t25 + o.BH_B - o.BH_C * t25
    check("pK_bh(25 degC) reproduces Bates & Allen 7.762", abs(pk25 - 7.762) < 5e-4, f"{pk25:.4f}")
    pkw25 = o.W_A / t25 - o.W_B + o.W_C * t25
    check("pK_w(25 degC) ~ 13.995", abs(pkw25 - 13.995) < 2e-3, f"{pkw25:.4f}")
    eps25 = o.MM_0 - o.MM_1 * 25 + o.MM_2 * 625 - o.MM_3 * 15625
    check("eps_r(25 degC) reproduces Malmberg & Maryott 78.30", abs(eps25 - 78.30) < 0.01, f"{eps25:.3f}")
    a25 = o.debye_hueckel_a(25.0)
    check("A(25 degC) ~ 0.511 (molar scale)", abs(a25 - 0.5115) < 1e-3, f"{a25:.5f}")

    # Oracle agrees with the helper re-implementation
    check("helper variant(all twists) == oracle", all(variant(*p) == o.compute(*p) for p in pts))

    # Twist integrity
    t1 = max(abs(o.compute(*p) - variant(*p, charge_aware=False)) for p in pts)
    t1b = max(abs(o.compute(*p) - variant(*p, activity=False)) for p in pts)
    t2 = max(abs(o.compute(*p) - variant(*p, thermal=False)) for p in pts)
    check("T1 charge-type activity twist changes output (> 100x tolerance)", t1 > 0.1 and t1b > 0.1,
          f"max diff vs acetate-type {t1:.4f}, vs ideal {t1b:.4f}")
    check("T2 thermal twist changes output (> 100x tolerance)", t2 > 0.1, f"max diff {t2:.4f}")
    check("T2 exactly neutral at value_d = 25",
          all(o.compute(p[0], p[1], p[2], 25.0) == variant(p[0], p[1], p[2], 25.0, thermal=False) for p in pts))
    up = all(o.compute(p[0], p[1], 0.2, p[3]) > o.compute(p[0], p[1], 0.0, p[3]) for p in pts)
    check("salt raises pH everywhere (charge-type sign)", up)
    s0 = o.compute(0.01, 0.01, 0.2, 0.0) - o.compute(0.01, 0.01, 0.0, 0.0)
    s50 = o.compute(0.01, 0.01, 0.2, 50.0) - o.compute(0.01, 0.01, 0.0, 50.0)
    check("salt effect grows with temperature (A(T) interaction)", s50 - s0 > 0.005,
          f"salt shift {s0:.4f} at 0 degC vs {s50:.4f} at 50 degC")
    cb = max(abs(o.compute(*p) - variant(*p, full_balance=False)) for p in pts)
    check("full charge balance matters somewhere in domain", cb > 0.01, f"max diff {cb:.4f}")

    # Golden data
    g = load_golden()
    cats = [c["category"] for c in g["cases"]]
    check(">=2 discriminating edge cases", cats.count("discriminating") >= 2)
    check(">=2 control edge cases", cats.count("control") >= 2)
    check(">=1 boundary edge case", cats.count("boundary") >= 1)
    check("golden inputs neutral and in range", all(
        sorted(c["input"]) == sorted(INPUT_KEYS)
        and all(INPUT_RANGES[k][0] <= c["input"][k] <= INPUT_RANGES[k][1] for k in INPUT_KEYS)
        for c in g["cases"]))
    check("golden expected outputs match oracle",
          all(o.oracle(c["input"]) == c["expected"] for c in g["cases"]), GOLDEN_PATH)

    print(json.dumps({"failures": failures}))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
