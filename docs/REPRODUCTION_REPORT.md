# Ponomarev 2016 — equation-to-code map, conditions, and results

**Implementation and reproducibility project: Kenshin Kotari.**

## 1. Source version and scope

**[P16] Anton Ponomarev**, *Nonlinear Predictor Feedback for Input-Affine Systems with Distributed Input Delays*, DOI [10.1109/TAC.2015.2496191](https://doi.org/10.1109/TAC.2015.2496191). Equation and page references in this report are to the inspected six-page **[arXiv:1601.00098v1](https://arxiv.org/abs/1601.00098v1), 1 January 2016**.

[P16, Section VI-A/B, p. 4] are analytical examples, without their own numerical figures. Our A/B runs are numerical illustrations with explicitly chosen parameters. [P16, VI-C, p. 5, Fig. 1] is the published pendulum simulation. These categories are not conflated.

We re-executed the existing physical-plant solvers from baseline `a491223df3ebc583e967fa81e4625d18afaf5824`, without editing the plant, feedback, history lookup, or test assertions. Our own implementation's limitations remain part of the record. A mismatch with earlier repository prose is not automatically a defect in the paper.

## 2. Equation-to-code map

| [P16] location | Source role | Executed function / use |
|---|---|---|
| **(4), (5)–(14)**, p. 1 | General input-affine plant and assumptions | Class being illustrated; no exhaustive assumption verification is claimed. |
| **(31)–(34)**, p. 2; **(46)–(47)**, p. 3 | Reconstruct predictor from current physical state and input history | Spatial loops in `run_example_a` and `run_example_c`. |
| **(39)–(42)**, p. 3 | Transformed dynamics and effective gain B | A/C beta loops; **(39) is not independently integrated to generate x**. |
| **(54)–(59), (64)**, p. 3 | Conditional stability and attraction-region bounds | Theorem context only; not established by finite trajectories or terminal test thresholds. |
| **(73)–(74)**, p. 4 | Scalar physical plant and feedback | [`python/src/example_a.py::run_example_a`](../python/src/example_a.py), MATLAB `run_example_a`. |
| **(75), (83), (85)–(88)**, p. 4 | Cascade physical plant, V, feedback, simplified transform | [`python/src/example_b.py::run_example_b`](../python/src/example_b.py), MATLAB `run_example_b`. |
| **(76)–(81)** versus **(86)–(88)**, p. 4 | Equivalent algebraic cascade expressions | `test_cross_validation.py::TestExampleBAnalyticalVsNumerical`; algebraic identity check, not independent plant simulation. |
| **(89)–(93), (104)**, p. 5 | Pendulum physical plant, Jacobian, gain dynamics, printed feedback | [`python/src/example_c.py::run_example_c`](../python/src/example_c.py), MATLAB `run_example_c`. |
| **(71), (98), (104)**, pp. 4–5 | General gradient feedback versus printed specialization | Normalization issue discussed below; the code uses **(104) literally**, not an unannounced replacement. |

## 3. Original equations and chosen conditions

### A. Scalar example

[P16, **(73)**]:

```math
\dot x=f(x)+b_0u(t)+b_1u(t-h)+\int_{-h}^0 b_{int}(\theta)u(t+\theta)\,d\theta.
```

[P16, **(74)**]: `u=-(f(y)+y)/B(y,u_t)`. The history-dependent predictor and gain are obtained using **(33)–(34), (40)–(42)**, rather than advancing the ideal `y_dot=-y` as the physical plant.

**Repository choices, not numerical values supplied in VI-A:** `f=sin`, `b0=1`, `b1=0.5`, constant `bint=1`, `h=0.5`, `x(0)=1`, constant negative-time input 0, `T=10`, and dt=0.002 or 0.001.

### B. Explicit cascade example

[P16, **(75)**]:

```math
\dot x_1=x_2^2+u(t-h),\qquad \dot x_2=x_2+u(t).
```

The code evaluates **(86)–(88)**: `z1=x1-x2+integral(u)`, `z2=x2`, `u=-2*x2+u_tilde`. With `q=z1+z2*(z2-2)/2`, **(83), (85)** give `V=q^2+z2^2` and `u_tilde=-(2*q*(z2-1)+2*z2)`. It then advances the physical right-hand side **(75)**, not cascade **(82)**.

**Repository choices, not a numeric simulation specification in VI-B:** `h=1`, `x(0)=[1,1]`, negative-time input 0, `T=15`, dt=0.002 or 0.001. The paper's remark after **(88)** explicitly says its proofs do not warrant global asymptotic stability of the original system for that example; we do not promote our finite run into that missing guarantee.

### C. Published pendulum example

[P16, **(89)**]:

```math
\dot x_1=x_2,\qquad \dot x_2=\sin(x_1)+u(t)+u(t-h),\qquad h=\pi/4.
```

From **(46)–(47)** with **(90)**, the predictor obeys `xi'=[xi2, sin(xi1)+u(t+s-h)]`, `xi(0)=x(t)`. From **(91)–(93)**, `beta'=A(xi)*beta`, `beta(0)=[0,1]`, `A=[[0,1],[cos(xi1),0]]`. The gain is `B=[0,1]+beta(h)` from **(40)**, and **(104)** is `u=-B.T@y`, `y=xi(h)`.

**From Fig. 1 and the paragraph after (104):** `x1(0)` in `{pi/2, pi, 3*pi/2}`, `x2(0)=pi/2`, initial input 1, and Euler step 0.01. The plotted horizon is 10. The additional central-angle runs at dt=0.02 and 0.005 are our refinement experiment.

**Source normalization observation, not a solver change:** inserting `v0=y.T@V@y` from **(98)** into **(71)** yields `-2*k*B.T@V@y`. With the stated `V=I,k=1` this differs by 2 from printed **(104)**. The existing solver uses **(104)** as written. Gain normalization may explain the presentation, but no undocumented reconciliation or proof correction is claimed.

## 4. Numerical contract

Forward Euler advances both physical states and A/C predictor/gain spatial equations. The spatial count is `N_h=round(h/dt)` and `ds=h/N_h`. Input lookup rounds `t_query/dt` to the nearest grid index and clamps to array bounds. At h=pi/4, dt=0.01 the physical delayed-input grid lag is 79 steps (0.79), while the predictor spatial interval still totals pi/4. This is a declared approximation, not exact delayed timestamp evaluation.

A/B history integrals use trapezoidal sums. Before the new command is stored, their current endpoint slot is still preallocated zero; the physical A integral is then recomputed after storing the command. A estimates f' using central differences with epsilon=1e-7 and sets `u=0` when `abs(B)<1e-10`. That guard is retained; its activation is not logged by the old solver, so this audit does not claim it was inactive. Final control values are held copies. These details matter when interpreting errors and comparing against ideal continuous-time formulas.

## 5. All measured runs

Data are **our local executions**, Python 3.13.5 / NumPy 2.3.5 / SciPy 1.17.0. See [summary.csv](../results/reporting-audit/summary.csv) and [provenance.json](../results/reporting-audit/provenance.json). All nine planned cases completed with finite arrays. No case was dropped for an unfavorable result. The runner and CI artifact retain the full run JSON and trajectories.

`norm(x(T))` is the Euclidean norm (absolute value in the scalar case). `max abs(u)` is the maximum over stored time nodes. Extra decimal places identify a recorded run, not certified continuous-solution accuracy.

| Case / source | dt | T | Final physical state | norm(x(T)) | u(0+) |
|---|---:|---:|---|---:|---:|
| A **(73)–(74)** | 0.002 | 10 | 0.0000578801874 | 0.0000578801874 | -1.110076161 |
| A **(73)–(74)** | 0.001 | 10 | 0.0000583660303 | 0.0000583660303 | -1.110228603 |
| B **(75), (85)–(88)** | 0.002 | 15 | [-0.000478021263, -0.0000939900977] | 0.000487173959 | -4.0 |
| B **(75), (85)–(88)** | 0.001 | 15 | [-0.000479043651, -0.0000941871341] | 0.000488215153 | -4.0 |
| C **(89), (104)**, x1(0)=pi/2 | 0.01 | 10 | [0.0236523098, -0.0127597473] | 0.0268745774 | -7.467661123 |
| C **(89), (104)**, x1(0)=pi | 0.01 | 10 | [0.0392582942, -0.0211772422] | 0.0446059329 | -6.811116128 |
| C **(89), (104)**, x1(0)=3pi/2 | 0.01 | 10 | [0.0571139653, -0.0308053963] | 0.0648920448 | -8.860879088 |
| C, central angle, refinement | 0.02 | 10 | [0.0391261633, -0.0210209478] | 0.0444155030 | -6.829414298 |
| C, central angle, refinement | 0.005 | 10 | [0.0393781337, -0.0212642890] | 0.0447527362 | -6.802301875 |

For B, stored final V from **(83)** is 2.3316478646e-8 (dt=0.002) and 2.3417036804e-8 (dt=0.001). Across these runs maximum absolute input equals the magnitude of the initial input shown above; the raw metadata contains that measurement.

**Interpretation:** A and B approach small residuals under the declared repository choices. C decreases substantially over the displayed horizon, but all three dt=0.01 terminal norms exceed 0.01. The old README's stricter threshold is contradicted by these executions. This is neither proof of asymptotic convergence nor proof of instability.

For central-angle C, the maximum componentwise physical-state difference relative to the linearly interpolated dt=0.005 run decreases from **0.0438882033** (dt=0.02) to **0.0309089373** (dt=0.01). See [refinement-comparison.json](../results/reporting-audit/refinement-comparison.json). This is agreement with a finer run of the same solver, not an independently known exact solution. A terminal norm need not decrease when the step decreases. No pixel-based fit to Fig. 1 has been performed.

## 6. What the tests actually check

[CI run 34175280499](https://github.com/KK1182112KK/ponomarev-2016-reproduction/actions/runs/34175280499) executed **32 existing Python tests: 32 passed**. The [raw log](../results/reporting-audit/tests.log) and [test record](../results/reporting-audit/test-results.json) are committed; full JUnit XML is in the CI artifact. That run also completed all nine numerical cases. Source hashes match our local solver sources; terminal norms agreed to approximately floating-point roundoff across the two recorded environments. This does not certify arbitrary platforms.

The pendulum terminal tests in `test_example_c.py` use **0.1**, not 0.01. `TestEulerVsRK4ExampleC` replays the **same recorded Euler command sequence** through an RK4 physical-plant step; it is not an independently closed-loop recomputation of the controller. The B cascade identity compares algebraically equivalent formulas using the same integral. None is a MATLAB/Python trajectory equivalence certificate. MATLAB was not executed in this reporting audit.

The CI source head was `ece0a41505184ccb7715bedf82ebde8b31a42fbc`; its test-merge checkout was `46affcfa6c402ff571751af001c79ab6bbbd6365`. CI recorded Python 3.13.15 / NumPy 2.5.3 / SciPy 1.18.1. These are execution records, not assertions that later documentation commits independently ran the tests.

## 7. Reproduction and remaining limits

```bash
pip install numpy scipy matplotlib pytest
python python/run_reporting_audit.py
python -m pytest python/tests -q
```

The runner regenerates all full CSV trajectories, logs, and source/environment metadata. Committed summaries retain all cases; full trajectories are supplied in the CI artifact with 30-day retention and can be regenerated. [METHODS.md](METHODS.md) specifies remaining numerical approximations. The original conditional theorem, its applicability, a certified attraction region, the original author's actual code, and a quantitatively exact figure reproduction have not been established by this project. The original theory remains credited to Ponomarev; these implementations, comparisons, and reports are the reproduction contribution.
