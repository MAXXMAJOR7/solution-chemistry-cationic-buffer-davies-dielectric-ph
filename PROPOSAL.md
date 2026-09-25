# Blackbox Task 48 - Cationic-Acid Buffer pH Under a Charge-Type Activity Inversion and Dielectric-Driven Temperature Laws

**Domain:** Chemistry - analytical / solution chemistry (acid-base equilibria, electrolyte activity theory, thermodynamics of ionisation)
**Tier:** Med (floor: 3h, 2 genuine twists, 3 substantive steps)
**Honest counts:** 4 substantive steps (plus a trivial rounding step), 2 genuine twists

**Repo name:** `solution-chemistry-cationic-buffer-davies-dielectric-ph`
**GitHub description:** Reverse-engineer the activity-scale pH of an amine/ammonium buffer in a salt background from 0 to 50 °C. A charge-symmetric acid makes added salt raise the pH, and the Debye-Hückel slope must be built from water's temperature-dependent permittivity. Requires electrolyte activity theory, exact speciation, and literature thermodynamics of ionisation.

---

## Task description (for reviewers)

The hidden function maps four real-valued inputs to the pH of an aqueous equilibrium system. Two literature effects combine inside a full speciation solve. The first is an electrolyte activity correction whose sign depends on the charge type of the acid. The second is a set of temperature laws, one of which is not tabulated anywhere and must be built from first principles plus a permittivity correlation. A solver must identify the system as a buffer, work out why inert salt pushes the pH *up*, and recover the exact temperature dependences to within 0.001 pH. Expected solver knowledge: graduate solution chemistry, including Debye-Hückel/Davies theory and the Harned-Robinson treatment of ionisation constants.

---

## 1. Domain & Algorithm

Solution chemistry. The quantity is the **activity-scale pH = −log₁₀ a_H+** of a buffer made from triethanolammonium chloride (BH⁺Cl⁻, "TEA·HCl") and triethanolamine (B) in water. The buffer contains an inert NaCl background and sits at a temperature between 0 and 50 °C. The oracle solves the complete equilibrium (buffer acid, water autoionisation, electroneutrality) with ionic strength made self-consistent. Single-ion activity coefficients come from the Davies equation, with the Debye-Hückel slope evaluated at temperature.

## 2. Core Method

- **Acid dissociation constant of BH⁺ vs T.** Bates & Allen (1960), *J. Res. NBS* 64A, 343, eq (5): −log K_bh = 1341.16/T + 4.6252 − 0.0045666 T, valid 0-50 °C (Harned-Robinson form).
- **Ion product of water vs T.** Harned & Robinson (1940), *Trans. Faraday Soc.* 36, 973: log K_w = −4470.99/T + 6.0875 − 0.01706 T.
- **Relative permittivity of water.** Malmberg & Maryott (1956), *J. Res. NBS* 56, 1: ε_r = 87.740 − 0.40008 t + 9.398×10⁻⁴ t² − 1.410×10⁻⁶ t³ (t in °C).
- **Debye-Hückel limiting slope from first principles.** ln γ = −z²e²κ/(8πε₀ε_r kT), with κ² = 2N_A e² I/(ε₀ε_r kT) (Debye & Hückel 1923; Robinson & Stokes, *Electrolyte Solutions*, ch. 9). A = slope/ln 10, in (mol L⁻¹)^−½, from exact SI constants and CODATA 2022 ε₀.
- **Davies equation** (Davies 1962, *Ion Association*): log₁₀ γ = −A z² (√I/(1+√I) − 0.3 I).
- **Speciation.** Electroneutrality h + [BH⁺] + [Na⁺] = [Cl⁻] + [OH⁻], with K_bh = a_H a_B / a_BH = h[B]/[BH⁺] (γ_H and γ_BH cancel; B is neutral) and K_w = γ²h[OH⁻].

## 3. Twists (2 genuine)

**T1 - Charge-type activity inversion (the salt raises the pH).** The textbook buffer, and the one every student computes with activities, is acetate: HA ⇌ H⁺ + A⁻. There K_a' = K_a/γ², and adding inert salt *lowers* the pH (pH = pK + log(b/a) + log γ). For the cationic acid BH⁺ ⇌ B + H⁺, the charges are the same on both sides. The activity coefficients cancel in the concentration quotient, but the operational pH still carries −log γ_H. The oracle therefore gives pH = pK + log(b/a) − log γ + (speciation terms), and salt *raises* the pH. The Davies turnover term (−0.3 I) makes the effect non-monotone in √I. A solver who applies the familiar acetate-type correction is off by up to 0.28 pH. Ignoring activity costs up to 0.14 pH.

**T2 - Temperature enters through two independent literature laws, one built from the solvent permittivity.**
(a) The buffer pK moves by ~−0.02 K⁻¹ via the Harned-Robinson three-term law, which curves (ΔC_p° ≠ 0). A constant-ΔH van't Hoff extrapolation from 25 °C misses by up to 0.034.
(b) The Debye-Hückel slope A(T) ∝ (ε_r T)^−3/2 is recomputed from Malmberg-Maryott ε_r(t) and CODATA constants. It rises from 0.492 at 0 °C to 0.537 at 50 °C, so the salt effect of T1 grows with temperature (an interaction). Freezing A at its 25 °C value fails 37/51 golden cases.
Water's K_w(T) is a third temperature law. Here it acts only through the charge balance.

Both twists are I/O-isolable:
- value_d = 25 makes T2 exactly neutral (verified bit-for-bit against the frozen variant).
- value_c = 0 with small value_a, value_b drives T1 towards neutrality (control cases). value_c alone switches it on.
- The T1 × T2 interaction appears only when value_c > 0 and value_d ≠ 25.

## 4. Multi-Step Pipeline

| Step | Computation | Inputs involved |
|---|---|---|
| S1 | pK_bh(T), pK_w(T), ε_r(t) → A(T) (T2) | value_d |
| S2 | Davies γ(I, A) for monovalent ions (T1 magnitude) | via I (value_a, value_c) and A |
| S3 | Solve electroneutrality for h by bisection on log h, iterating I = [Cl⁻] + [OH⁻] to self-consistency | value_a, value_b, value_c |
| S4 | pH = −log₁₀(γ h), where the charge type sets the sign of the salt effect (T1) | all |
| (S5) | round to 4 decimals | trivial, not counted |

The steps are separable and each is substantive. The solve in S3 is a real step: at 200:1 base excess and 0 °C, [OH⁻] consumes ~8% of the BH⁺. Dropping it moves the pH by up to 0.07 across the domain and by 0.036 at the golden boundary case.

## 5. Input Schema

```json
{
  "value_a": "float [0.0005, 0.1]",   // triethanolammonium chloride, mol L^-1
  "value_b": "float [0.0005, 0.1]",   // triethanolamine (free base), mol L^-1
  "value_c": "float [0, 0.2]",        // NaCl, mol L^-1
  "value_d": "float [0, 50]"          // temperature, degC (validity range of Bates & Allen eq 5)
}
```

## 6. Output Schema

```json
{ "output": "float, 4 decimals" }     // pH (activity scale); range ~5.1 to 10.6
```

## 7. Hardness & Twists

These baselines are measured on the 51 golden cases (`.scratch/baselines.py`, tolerance 0.001 pH):

| Model | Pass | Median abs err | Max abs err |
|---|---|---|---|
| Henderson-Hasselbalch, 25 °C pK | 0/51 | 0.267 | 0.649 |
| Acetate-type activity correction (salt lowers pH) | 0/51 | 0.223 | 0.284 |
| Ideal solution (γ = 1), exact T laws | 0/51 | 0.112 | 0.142 |
| All T laws frozen at 25 °C (T2 off) | 4/51 | 0.259 | 0.570 |
| Constant-ΔH van't Hoff pK(T) | 18/51 | 0.0021 | 0.034 |
| A frozen at 25 °C (only A(T) off) | 14/51 | 0.0023 | 0.0068 |
| No h / [OH⁻] in charge balance | 35/51 | 0.0004 | 0.036 |
| Full model (example solver, 12 probes) | 51/51 | 4e-6 | - |

Why this is hard:
- **The sign of the salt effect is a trap.** Almost all worked examples of activity-corrected buffers use anionic (HA/A⁻) or phosphate acids, where salt lowers the pH. A solver who hypothesises "buffer + Davies" without reasoning about charge type fits a curve with the wrong sign.
- **A(T) is not a lookup.** Recovering the small but tolerance-dominant T × salt interaction requires the (ε_r T)^−3/2 scaling and an accurate ε_r(t). A linear A(T) or a 25 °C constant fails.
- **The pK(T) law is curved.** A van't Hoff line from 25 °C fails at both ends. The three-parameter Harned-Robinson form must be identified.
- **Linear or additive surrogates fail.** The output is logarithmic in value_b/value_a but bends at extreme ratios (water term, free h). It is non-linear in √I, and the salt effect is multiplied by a temperature factor.
- **Cold read:** four neutral floats and an output in ~5-11 suggest "some log quantity". Recognising a buffer is plausible from the +1/decade ratio law. It does not reveal the acid's identity, its charge type, the temperature laws, or the activity model.

## 8. Constants & Citations

Every numeric literal in `oracle/implement.py` is whitelisted with its justification in `scripts/verify_task.py` (AST scan).

| Constant | Value | Status | Source |
|---|---|---|---|
| e | 1.602176634×10⁻¹⁹ C | EXACT | SI 2019 defining constant |
| k_B | 1.380649×10⁻²³ J K⁻¹ | EXACT | SI 2019 defining constant |
| N_A | 6.02214076×10²³ mol⁻¹ | EXACT | SI 2019 defining constant |
| ε₀ | 8.8541878188×10⁻¹² F m⁻¹ | CITED | CODATA 2022 (NIST) |
| 273.15, 1000 | - | EXACT | °C → K; L per m³ |
| 1341.16, 4.6252, 0.0045666 | - | CITED | Bates & Allen 1960, eq (5) |
| 4470.99, 6.0875, 0.01706 | - | CITED | Harned & Robinson 1940 |
| 87.740, 0.40008, 9.398e-4, 1.410e-6 | - | CITED | Malmberg & Maryott 1956 |
| 0.3 | - | CITED | Davies 1962 |
| 8, 2 (in DH slope) | - | EXACT | Debye-Hückel theory (8πε kT; κ² = 2N_A e² I/εkT) |
| −14, −1 | - | STRUCTURAL | bisection bracket on log₁₀h (pH 1 to 14 contains the domain's 5.1 to 10.6) |
| 200, 60 | - | STRUCTURAL | fixed iteration counts (bisection converges to machine precision; fixed point contracts because I depends on h only through the small [OH⁻] term) |
| 50, 0 | - | STRUCTURAL | input-domain check (Bates & Allen validity 0-50 °C) |

### Citations fetched & verified (quote < 15 words, location)

| Item | Verbatim quote | Where fetched |
|---|---|---|
| pK_bh(T) | "−log K_bh = 1341.16/T + 4.6252 − 0.0045666 T (5)" | Bates & Allen 1960, §4 eq (5), Europe PMC full text PMC5287088 |
| pK_bh(25 °C) | "The value of −log K_bh at 25° C (7.762)" | same, §3 (oracle gives 7.7619) |
| Validity | "Between 0° and 50° C, K_bh is given by the expression" | same, §4 |
| K_w(T) | "log κ2 = log Kw = 6.0875 − 4470.99/T − 0.01706T" | arXiv:1405.1338, eq (2.18), citing Harned & Robinson form |
| ε_r(t) | "fit the equation ε = 87.740 − 0.40008t + 9.398(10⁻⁴)t² − 1.410(10⁻⁶)t³" | Malmberg & Maryott 1956, archive.org full text (jresv56n1p1)* |
| Davies | "−log f± = 0.5 z₁z₂ (√I/(1+√I) − 0.30 I)" | OECD-NEA TDB-2 guidelines (Grenthe), citing Davies 1962 |
| ε₀ | "8.854 187 8188 x 10⁻¹² F m⁻¹" | physics.nist.gov CODATA 2022 value page |

\*The archive.org OCR reads the linear coefficient as "0.4008". The paper's own 25 °C value, 78.30, is reproduced only with 0.40008 (78.303 vs 78.285), and 0.40008 is the value in all secondary reproductions, so 0.40008 is used.

**Citation-code match:** `BH_A, BH_B, BH_C = 1341.16, 4.6252, 0.0045666`; `W_A, W_B, W_C = 4470.99, 6.0875, 0.01706`; `MM_0..MM_3 = 87.740, 0.40008, 9.398e-4, 1.410e-6`; `DAVIES_B = 0.3`; `EPS_0 = 8.8541878188e-12`. `verify_task.py` asserts the literature anchors: pK_bh(25) = 7.762, ε_r(25) = 78.30, A(25) ≈ 0.511.

**Honest notes (CONCERN-level, documented):**
- **Davies beyond its comfort zone.** NEA quotes Davies as reliable to ~0.1 mol kg⁻¹. The domain reaches I = 0.3. The oracle is a deterministic model, not a claim of experimental accuracy there.
- **Molar vs molal.** Concentrations are treated as mol L⁻¹ throughout, and A is computed on the molar scale (the κ² expression with 1000 L m⁻³). Density corrections are deliberately omitted.
- **Single-ion activity.** pH = −log(γh) uses a Davies single-ion coefficient, which is the conventional operational choice (as in PHREEQC and MINEQL). It is a convention, not a measurable quantity.
- **Secondary source for K_w.** The Harned & Robinson (1940) primary text was not open-access. The coefficients were verified in an open paper quoting that law. They reproduce pK_w(25 °C) = 13.995.
- **Sanity check against Bates & Allen's own data.** Table 1 reports pwH (which includes γ_Cl). For m₁ = 0.09909, m₂ = 0.10481 at 25 °C, Davies predicts 8.000 vs the measured 8.047. The sign and magnitude of the salt effect are right. The residual reflects Davies vs their two-parameter Hückel fit.

## 9. Edge Cases

| id | value_a | value_b | value_c | value_d | output | Rationale |
|---|---|---|---|---|---|---|
| edge_control_1 | 0.001 | 0.001 | 0 | 25 | 7.7769 | T2 neutral; T1 small (pK + 0.015) |
| edge_control_2 | 0.002 | 0.02 | 0 | 25 | 8.7820 | T2 neutral; textbook pK + 1 (+0.02 activity) |
| edge_control_3 | 0.0005 | 0.0005 | 0 | 25 | 7.7720 | Approaches pK(25 °C) = 7.7619 |
| edge_discrim_1 | 0.01 | 0.01 | 0.2 | 25 | 7.8904 | T1: salt raises pH by 0.13 (acetate-type model: 7.6334) |
| edge_discrim_2 | 0.001 | 0.001 | 0 | 0 | 8.3025 | T2 alone: +0.53 |
| edge_discrim_3 | 0.001 | 0.001 | 0 | 50 | 7.3151 | T2 alone: −0.46 |
| edge_discrim_4 | 0.01 | 0.01 | 0.2 | 50 | 7.4346 | T1 × T2: larger A(T) → larger activity shift (+0.135 over ideal) |
| edge_discrim_5 | 0.01 | 0.01 | 0.2 | 0 | 8.4113 | T1 × T2: smaller A(T) (+0.124 over ideal) |
| edge_boundary_1 | 0.0005 | 0.1 | 0 | 0 | 10.5639 | Max base excess, 0 °C: water term shifts pH by −0.036 |
| edge_boundary_2 | 0.1 | 0.0005 | 0 | 50 | 5.1202 | Max acid excess, 50 °C: free h matters (+0.009) |
| edge_boundary_3 | 0.1 | 0.1 | 0.2 | 50 | 7.4416 | Max ionic strength (I = 0.3), largest activity shift |

The golden file adds 40 seeded random cases (seed 48), for 51 cases in total.

---

## GATE 1: Screening

| Test | Result | Notes |
|---|---|---|
| Algebraic collapse | PASS | In the mid range the output ≈ pK(value_d) + log(value_b/value_a) − log γ(value_a + value_c, A(value_d)). The ratio term is additive, but γ couples value_a, value_c and value_d non-separably. The exact charge balance adds a further non-separable term at the extremes. |
| Domain recall | PASS (CONCERN) | An expert may say "buffer with activity correction". The specific acid, its charge type, the Harned-Robinson laws, and A(T) from ε_r(t) are not named by the schema. No single published method is this composition. |
| Twist survives | PASS | T1 changes the sign of an S-shaped √I function. T2 is a curved 1/T + T law times an (ε_r T)^−3/2 factor. Neither simplifies away. |
| Genuine twist | PASS (2) | T1 contradicts the standard activity-corrected buffer result. T2 replaces 25 °C constants with two independent literature laws that interact with T1. K_w(T) and the charge balance are pipeline rigour and are not counted. |
| Two-trap | PASS | value_c and value_d are independently probeable. value_d = 25 isolates T1, and value_c = 0 (dilute) isolates T2. |
| Tier fit | PASS | 4 steps ≥ 3; 2 twists ≥ 2 (Med). |
| PhD authenticity | PASS | The Bates-school pH-standard work rests on exactly these distinctions (charge type, T laws). Misapplying the acetate-type correction to amine buffers (Tris, TEA, HEPES-type) is a real laboratory error. |

**Verdict: PASS** (one documented CONCERN on domain recall of the buffer family)

## GATE 2: Verification

| Check | Result |
|---|---|
| Twist integrity (freeze-one tests) | PASS: over 300 random points, T1 max effect 0.28 (vs acetate-type) / 0.14 (vs ideal), T2 max effect 0.52 (`verify_task.py`) |
| I/O isolation | PASS: value_d = 25 ⇒ T2 bit-exact neutral; salt raises pH at every sampled point; salt effect 0.080 at 0 °C vs 0.088 at 50 °C |
| No invented constant | PASS: AST scan finds only whitelisted literals |
| Recall-reconstructability | PASS: the example solver recovers the function with 12 probes (≪ 50) once the family is hypothesised |
| Citations fetched & verified | PASS: all seven sources fetched (table above), one OCR discrepancy resolved by an internal anchor |
| Citation-code match | PASS |
| Output/schema neutralisation | PASS: key "output", float, value_a to value_d |
| Tier floor | PASS: 4 steps / 2 twists vs Med 3 / 2 |
| Edge cases ≥ 5 | PASS: 11 total (3 control, 5 discriminating, 3 boundary) |
| PhD-level QA | PASS: needs activity theory with charge-type reasoning, dielectric-based DH slope, and Harned-Robinson thermodynamics together |

**Verdict: PASS**

---

## Files

```
task-48/
├── oracle/implement.py          reference oracle (hidden)
├── solution/example_solver.py   probe + fit solver (12 probes, 51/51)
├── grader/grade_submission.py   tolerance 0.001 pH absolute, probe budget 200
├── golden/test_data.json        51 cases (11 edge + 40 random)
├── scripts/helpers.py           shared utilities, twist-freezing variant()
├── scripts/build_golden.py      regenerates golden data
├── scripts/verify_task.py       automated GATE 2 checks
├── scripts/run_tests.py         verify + oracle self-grade + solver grade
└── PROPOSAL.md
```

## References

- Bates, R. G. & Allen, G. F. (1960). Acid dissociation constant and related thermodynamic quantities for triethanolammonium ion in water from 0 to 50 °C. *J. Res. Natl. Bur. Stand.* 64A, 343-346. doi:10.6028/jres.064A.033 (PMC5287088)
- Harned, H. S. & Robinson, R. A. (1940). A note on the temperature variation of the ionisation constants of weak electrolytes. *Trans. Faraday Soc.* 36, 973-978. doi:10.1039/TF9403600973
- Malmberg, C. G. & Maryott, A. A. (1956). Dielectric constant of water from 0° to 100 °C. *J. Res. Natl. Bur. Stand.* 56, 1-8. doi:10.6028/jres.056.001
- Davies, C. W. (1962). *Ion Association*. Butterworths, London.
- Debye, P. & Hückel, E. (1923). Zur Theorie der Elektrolyte. *Phys. Z.* 24, 185-206.
- Robinson, R. A. & Stokes, R. H. (1959). *Electrolyte Solutions*, 2nd ed. Butterworths (ch. 9; app. 7.1).
- Grenthe, I. et al. TDB-2: Guidelines for the extrapolation to zero ionic strength. OECD Nuclear Energy Agency.
- CODATA 2022 recommended values, NIST. https://physics.nist.gov/cuu/Constants/
- Open-access verification sources: Europe PMC PMC5287088; arXiv:1405.1338; archive.org jresv56n1p1.
