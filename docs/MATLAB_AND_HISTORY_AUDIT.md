# MATLAB execution, cross-language comparison, and history audit

Implementation/reproducibility project: **Kenshin Kotari**. Investigation date: 2026-09-08 UTC. Original theory: Anton Ponomarev. Equation numbers below follow **arXiv:1601.00098v1**, not presumed journal renumbering. Development and investigation were AI-assisted.

## Executed evidence, not an inferred pass

The previously unexecuted MATLAB implementation has now been run on GitHub Actions with **MATLAB 26.1.0.3346908 (R2026a) Update 5**. All **18 discovered existing tests passed**; all **six declared simulation cases completed**. No original solver or test assertion was changed.

- Audited production baseline: `bc3d17fe37015236997f6535cc3b6a885896bd8d`.
- Audit harness head: `c498ddb74fcaf9a9672b6723fd8a75ac7ac2c621`.
- Actual PR test-merge checkout: `19b92a163fa2928b2a8c8453a91c4974e8f63f58`.
- [MATLAB run 34180001445](https://github.com/KK1182112KK/ponomarev-2016-reproduction/actions/runs/34180001445), **attempt 2**; artifact `10038656300`, `matlab-execution-audit`.
- ZIP SHA256: `32ae89c2ada528f7a2951d0a21296efa036718bc9949a7c8d4a409c0f1445c05`.
- Attempt 1 failed during MATLAB installation with `double free or corruption (fasttop)`. Tests were not executed in that attempt. A retry, without solver changes, produced the evidence above.

The artifact contains `audit.json`, `tests.csv`, `tests.xml`, `matlab.log`, `source-sha256.txt`, `checkout.txt`, and six full trajectory CSVs. Its retention is 90 days. The compact measured record and raw test table are also committed under [`results/matlab-history-audit/`](../results/matlab-history-audit/). Local Python 3.13.5 reruns matched all **21** source/test hashes in the MATLAB manifest before comparison.

## 1. What was compared

Example A implements physical Eq. **(73)** and feedback **(74)** through predictor **(31)-(34)** and gain **(40)-(42)**. Example B advances physical Eq. **(75)**, reconstructs **(86)-(87)**, and applies **(85), (88)**. Example C advances physical Eq. **(89)** and uses predictor **(46)-(47)**, gain **(40), (92)-(93)**, and feedback **(104)**. Neither language generates physical X from a pre-solved stable target.

These remain Forward Euler / nearest-grid-history implementations. Cross-language agreement tests the implementation equivalence for stated conditions, **not exact continuous-time reproduction**. Define the reported errors as the maximum absolute component difference over corresponding output nodes, separately for X and U. These are not Euclidean trajectory norms; CSV serialization also contributes rounding.

| Case | Conditions | MATLAB final physical norm | max componentwise X difference | max componentwise U difference |
|---|---|---:|---:|---:|
| A, scalar | h=.5, dt=.001, T=10, x0=1, past U=0, f=sin, b0=1, b1=.5, bint=1 | 0.0000583660302896886 | 1.09024e-13 | 1.44456e-11 |
| B, cascade | h=1, dt=.001, T=15, x0=[1,1], past U=0 | 0.000488215153355191 | 5.10703e-15 | 4.88498e-15 |
| C, initial angle pi/2 | h=pi/4, dt=.01, T=10, x2(0)=pi/2, past U=1 | 0.0268745773814305 | 5.32907e-15 | 5.32907e-15 |
| C, initial angle pi | same settings | 0.0446059329334678 | 5.10703e-15 | 4.88498e-15 |
| C, initial angle 3pi/2 | same settings | 0.0648920447794677 | 5.32907e-15 | 4.88498e-15 |
| C, half-grid diagnostic | **audit-only** h=.45, dt=.02, T=.1, x0=[pi,pi/2], past U=1 | 3.49975292387723 | **0.00178393518** | **0.0341838176** |

The standard cases agree closely. The half-grid case does not. This is retained, not discarded. A/B settings are repository illustrations, not published simulation parameters.

## 2. Half-grid rounding changes the discretization

The unchanged solvers choose `N_h=round(h/dt)`. For h=.45 and dt=.02, h/dt=22.5. Python rounds this to **22**, whereas MATLAB's default rounds it to **23**. Therefore predictor subdivisions and subsequent history queries differ. These are language rounding conventions, not random solver noise. The artifact records MATLAB's tie examples; the standard-library definitions are [Python round](https://docs.python.org/3/library/functions.html#round) and [MATLAB round](https://www.mathworks.com/help/matlab/ref/double.round.html).

The half-grid comparison is an additional numerical portability test, not a reproduction of Fig. 1. It rules out an unrestricted claim of cross-language identity for arbitrary parameter choices. It does not negate the agreement of the standard cases.

## 3. Strict negative-time history is not always preserved

[`run_history_probe.py`](../python/run_history_probe.py) calls each unchanged `_lookup_u` in `example_a.py`, `example_b.py`, and `example_c.py` with a deliberately discontinuous test history:

- dt=.01, requested time q=-.004;
- all negative-time stored inputs equal 1;
- the command at time zero equals 9.

All three helpers return **9**, whereas the strict negative-time history is **1**. Rounding q/dt=-.4 to zero crosses the startup boundary. This is a concrete failure to preserve the specified prehistory at a non-grid query. It does **not** establish that every default trajectory makes this particular query, nor that the algorithm accesses measured future physical states.

## 4. Eq. (86): the implemented quadrature is not the actual held-input integral

For physical Eq. **(75)** the transformation **(86)** is

```math
z_1(t)=x_1(t)-x_2(t)+\int_{t-h}^{t}u(s)\,ds.
```

The existing B code evaluates a trapezoidal sum before it computes/stores the new U[k]. Its upper endpoint is therefore an **unissued, preallocated zero slot**, with nonzero quadrature weight. A single endpoint has zero weight in the exact integral, but not in this numerical quadrature. The same pattern appears in the scalar predictor's initial integral.

An independent postprocessing probe uses the exact integral of the issued piecewise-constant history: for h=1, dt=.001, past U=0, it sums the preceding 1000 held bins. With x0=[1,1] and T=2, the maximum discrepancy in reconstructed z1 is **0.0020000000000000018**. At t=1:

| Quantity | Value |
|---|---:|
| Stored z1 | -0.09999366891399974 |
| z1 reconstructed from the actual ZOH input | -0.10199366891399975 |

This diagnoses the difference between the legacy quadrature convention and a strict ZOH-history interpretation of **(86)**. It is not a theorem counterexample or a claim that all finite-step computations are unusable. Algebraic agreement of two transforms that share a quadrature value cannot independently validate that value.

## 5. Interpretation and remaining work

The original physical RHS is present, the standard MATLAB/Python outputs agree, and the existing MATLAB tests pass. However, those tests do not remove the demonstrated non-grid history and quadrature-endpoint issues. These findings concern **this repository's numerical implementation**, not misconduct or invalidity of Ponomarev's theory.

No production solver or existing assertion is changed in this investigation. A causal timestamp lookup, an explicit endpoint convention, and regression tests for off-grid queries are appropriate next implementation changes; any such change requires its own before/after and refinement record. No certified region of attraction or asymptotic stability conclusion follows from the finite tests.

## Re-run

```bash
python python/run_history_probe.py
# Download/extract the cited MATLAB artifact first; this command checks source hashes.
python python/compare_matlab_audit.py --matlab-dir /path/to/extracted/artifact
```

The history probe is a **diagnostic recorder**, not a test that reports success when a defect exists; inspect its boolean fields. The comparison script reproduces all artifact cases, including the mismatching half-grid case. MATLAB execution is reproducible through `addpath('matlab'); run_cross_language_audit` and the new audit workflow. Full logs and trajectories remain in the artifact; compact evidence is committed separately. Later documentation commits are not retroactively described as the tested checkout.
