"""Regenerate golden/test_data.json from the oracle (deterministic, seeded)."""
import json
import math
import random

from helpers import GOLDEN_PATH, load_oracle

# (id, category, label, value_a, value_b, value_c, value_d)
EDGE_CASES = [
    ("edge_control_1", "control", "25 degC, dilute equimolar, no salt: thermal twist neutral, activity term small (pH ~ pK + 0.01)",
     0.001, 0.001, 0.0, 25.0),
    ("edge_control_2", "control", "25 degC, 1:10 buffer, no salt: thermal twist neutral, near textbook pK + 1",
     0.002, 0.02, 0.0, 25.0),
    ("edge_control_3", "control", "Same as control_1 at 0.0005 M: activity term vanishes further, output approaches pK(25 degC) = 7.7619",
     0.0005, 0.0005, 0.0, 25.0),
    ("edge_discrim_1", "discriminating", "25 degC, 0.2 M NaCl: charge-type activity twist RAISES pH by ~0.13 over ideal (acetate-type intuition predicts a ~0.13 fall)",
     0.01, 0.01, 0.2, 25.0),
    ("edge_discrim_2", "discriminating", "0 degC, dilute: thermal twist alone lifts pH by ~0.53",
     0.001, 0.001, 0.0, 0.0),
    ("edge_discrim_3", "discriminating", "50 degC, dilute: thermal twist alone lowers pH by ~0.46",
     0.001, 0.001, 0.0, 50.0),
    ("edge_discrim_4", "discriminating", "50 degC, high salt: A(T) from eps_r(T) enlarges the salt effect vs discrim_1",
     0.01, 0.01, 0.2, 50.0),
    ("edge_discrim_5", "discriminating", "0 degC, high salt: both twists fire, smaller A(T)",
     0.01, 0.01, 0.2, 0.0),
    ("edge_boundary_1", "boundary", "Corner: 200:1 base excess at 0 degC -> pH > 10.5, [OH-] consumes ~8% of BH+ (full charge balance needed)",
     0.0005, 0.1, 0.0, 0.0),
    ("edge_boundary_2", "boundary", "Corner: 200:1 acid excess at 50 degC -> lowest pH ~5.1, free h (~8e-6 M) no longer negligible",
     0.1, 0.0005, 0.0, 50.0),
    ("edge_boundary_3", "boundary", "Corner: maximum ionic strength (I = 0.3), equimolar 0.1 M, 50 degC: largest activity shift (~0.14)",
     0.1, 0.1, 0.2, 50.0),
]

N_RANDOM = 40
SEED = 48


def main():
    o = load_oracle()
    cases = []
    for cid, cat, label, a, b, c, d in EDGE_CASES:
        inp = {"value_a": a, "value_b": b, "value_c": c, "value_d": d}
        cases.append({"id": cid, "category": cat, "label": label, "input": inp, "expected": o.oracle(inp)})
    rng = random.Random(SEED)
    lo, hi = math.log10(0.0005), math.log10(0.1)
    for i in range(N_RANDOM):
        inp = {
            "value_a": round(10 ** rng.uniform(lo, hi), 6),
            "value_b": round(10 ** rng.uniform(lo, hi), 6),
            "value_c": round(rng.uniform(0.0, 0.2), 4),
            "value_d": round(rng.uniform(0.0, 50.0), 2),
        }
        cases.append({"id": f"rand_{i:02d}", "category": "random", "label": "log-uniform value_a, value_b sample",
                      "input": inp, "expected": o.oracle(inp)})
    doc = {
        "schema": {
            "input": {
                "value_a": "float [0.0005, 0.1]",
                "value_b": "float [0.0005, 0.1]",
                "value_c": "float [0, 0.2]",
                "value_d": "float [0, 50]",
            },
            "output": {"output": "float, 4 decimals"},
        },
        "seed": SEED,
        "cases": cases,
    }
    assert all(math.isfinite(c["expected"]["output"]) for c in cases)
    with open(GOLDEN_PATH, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    print(f"wrote {len(cases)} cases to {GOLDEN_PATH}")


if __name__ == "__main__":
    main()
