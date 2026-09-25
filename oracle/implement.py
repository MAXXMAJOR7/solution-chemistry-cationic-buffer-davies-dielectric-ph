"""Reference oracle (hidden from solver).

pH (= -log10 a_H+) of an aqueous triethanolamine / triethanolammonium chloride
buffer carrying an inert 1:1 background electrolyte (NaCl), at temperature t.

Pipeline (see PROPOSAL.md sections 2-4):
  S1  temperature laws        pK_bh(T) = 1341.16/T + 4.6252 - 0.0045666 T   (Bates & Allen 1960)
                              pK_w(T)  = 4470.99/T - 6.0875 + 0.01706 T     (Harned & Robinson 1940)
                              eps_r(t) = 87.740 - 0.40008 t + 9.398e-4 t^2 - 1.410e-6 t^3
                                                                          (Malmberg & Maryott 1956)
                              A(T)     = Debye-Hueckel slope from CODATA constants and eps_r(t)  (T2)
  S2  single-ion activity     log10 gamma = -A (sqrt I / (1 + sqrt I) - 0.3 I)   (Davies 1962)
  S3  speciation              exact charge balance
                                h + [BH+] + [Na+] = [Cl-] + [OH-],
                                [BH+] = C_T h / (h + K_bh)   (charge-symmetric acid: gammas cancel),
                                [OH-] = K_w / (gamma^2 h),
                              solved together with I = [Cl-] + [OH-]  (self-consistent)
  S4  operational pH          pH = -log10(gamma h)   -> the BH+ charge type makes salt RAISE pH  (T1)

Inputs (neutral names):
    value_a : triethanolammonium chloride, mol L^-1   [0.0005, 0.1]
    value_b : triethanolamine (free base), mol L^-1   [0.0005, 0.1]
    value_c : NaCl, mol L^-1                          [0, 0.2]
    value_d : temperature, degC                       [0, 50]
Output:
    output  : pH, rounded to 4 decimals
"""
import math

# --- SI exact defining constants and CODATA 2022 --------------------------------
E_CHARGE = 1.602176634e-19      # C, elementary charge (exact, SI 2019)
K_B = 1.380649e-23              # J K^-1, Boltzmann constant (exact, SI 2019)
N_A = 6.02214076e23             # mol^-1, Avogadro constant (exact, SI 2019)
EPS_0 = 8.8541878188e-12        # F m^-1, vacuum permittivity (CODATA 2022)
T_ZERO = 273.15                 # K, 0 degC (exact)
L_PER_M3 = 1000.0               # L m^-3 (exact): mol L^-1 -> mol m^-3

# --- Bates & Allen (1960), J. Res. NBS 64A, 343, eq (5): -log K_bh, 0-50 degC ----
BH_A, BH_B, BH_C = 1341.16, 4.6252, 0.0045666

# --- Harned & Robinson (1940), Trans. Faraday Soc. 36, 973: log K_w -------------
W_A, W_B, W_C = 4470.99, 6.0875, 0.01706

# --- Malmberg & Maryott (1956), J. Res. NBS 56, 1: eps_r(t), 0-100 degC ----------
MM_0, MM_1, MM_2, MM_3 = 87.740, 0.40008, 9.398e-4, 1.410e-6

# --- Davies (1962), Ion Association: linear-term coefficient ----------------------
DAVIES_B = 0.3

N_BISECT = 200                  # fixed iteration counts -> bit-for-bit determinism
N_FIXED_POINT = 60


def debye_hueckel_a(t_c):
    """Debye-Hueckel limiting slope A (log10 basis, (mol L^-1)^-1/2) at t_c degC."""
    temp = t_c + T_ZERO
    eps_r = MM_0 - MM_1 * t_c + MM_2 * t_c ** 2 - MM_3 * t_c ** 3
    ekt = EPS_0 * eps_r * K_B * temp
    # ln gamma = -z^2 e^2 kappa / (8 pi eps kT),  kappa^2 = 2 N_A e^2 I / (eps kT)
    slope_ln = E_CHARGE ** 2 / (8 * math.pi * ekt) * math.sqrt(2 * N_A * L_PER_M3 * E_CHARGE ** 2 / ekt)
    return slope_ln / math.log(10)


def davies_gamma(a_dh, ionic):
    s = math.sqrt(ionic)
    return 10 ** (-a_dh * (s / (1 + s) - DAVIES_B * ionic))


def compute(value_a, value_b, value_c, value_d):
    c_bh = float(value_a)
    c_b = float(value_b)
    c_na = float(value_c)
    t_c = float(value_d)
    if not (c_bh > 0 and c_b > 0 and c_na >= 0 and 0 <= t_c <= 50):
        raise ValueError("inputs out of domain")

    # S1: temperature laws
    temp = t_c + T_ZERO
    k_bh = 10 ** -(BH_A / temp + BH_B - BH_C * temp)
    k_w = 10 ** (-W_A / temp + W_B - W_C * temp)
    a_dh = debye_hueckel_a(t_c)

    c_tot = c_bh + c_b
    c_cl = c_bh + c_na
    ionic = c_cl
    h = gamma = None
    for _ in range(N_FIXED_POINT):
        # S2: Davies activity coefficient (all ions monovalent)
        gamma = davies_gamma(a_dh, ionic)

        # S3: charge balance, monotone increasing in h -> bisection on log10 h
        def residual(log_h):
            hh = 10 ** log_h
            return hh + c_tot * hh / (hh + k_bh) + c_na - c_cl - k_w / (gamma ** 2 * hh)

        lo, hi = -14.0, -1.0
        for _ in range(N_BISECT):
            mid = 0.5 * (lo + hi)
            if residual(mid) > 0:
                hi = mid
            else:
                lo = mid
        h = 10 ** (0.5 * (lo + hi))
        ionic = c_cl + k_w / (gamma ** 2 * h)

    # S4: operational pH on the activity scale
    return round(-math.log10(gamma * h), 4)


def oracle(inputs):
    return {"output": compute(inputs["value_a"], inputs["value_b"], inputs["value_c"], inputs["value_d"])}


if __name__ == "__main__":
    import json
    import sys
    print(json.dumps(oracle(json.loads(sys.stdin.read()))))
