# Ponomarev (2016): equation-linked numerical reproduction

**Reproduction implementation and verification project: Kenshin Kotari.** Original models, predictor design, and theorems are Anton Ponomarev's. Development and documentation were AI-assisted. Cite this software separately using [CITATION.cff](CITATION.cff).

**Source [P16]:** *Nonlinear Predictor Feedback for Input-Affine Systems with Distributed Input Delays*, DOI [10.1109/TAC.2015.2496191](https://doi.org/10.1109/TAC.2015.2496191). **All equation numbers here follow [arXiv:1601.00098v1](https://arxiv.org/abs/1601.00098v1), pages 1–6**, not an assumed renumbering of the journal version.

## Equations and measured outcomes

The default solvers advance the physical plant with Forward Euler and stored-input history; predictors are reconstructed at each control update. They do not generate the physical response from an independently stabilized target trajectory. They retain nearest-grid delay lookup and their original numerical conventions, not an exact continuous-time DDE implementation.

| Project example | Source equations | Actual recorded outcome |
|---|---|---|
| A — scalar distributed input | Plant **(73)**, feedback **(74)**, predictor/gain **(31)–(34), (40)–(42)** | Repository-chosen sine example: at dt=0.001, T=10, **abs(x)=0.0000583660303**. These numerical parameters are not specified in VI-A. |
| B — explicit cascade | Plant **(75)**, transforms **(76)–(88)**, feedback **(85), (88)** | At dt=0.001, T=15, **norm(x)=0.000488215153**. This is our numerical illustration of VI-B, not a published figure reproduction. |
| C — inverted pendulum | Plant **(89)**, predictor **(46)–(47)**, gain **(40), (92)–(93)**, feedback **(104)**, **Fig. 1** | At dt=0.01, T=10, the three initial angles give state norms **0.0268745774, 0.0446059329, 0.0648920448**. None is below 0.01 at T=10. |

[Full equation-to-code maps, all nine runs, conditions, and interpretation](docs/REPRODUCTION_REPORT.md) · [Measured CSV](results/reporting-audit/summary.csv)

The previous README's blanket claim that Example C reaches a norm below 1e-2 is **not supported by these runs**. The existing tests use a 0.1 terminal threshold for these cases; passing those tests does not establish the stronger README claim. A finite nonzero terminal norm is not evidence of instability. No exact pixel-level reproduction, global-stability proof, or MATLAB/Python equivalence is claimed.

## Execute the reproducible audit

[Open the audit notebook in Colab](https://colab.research.google.com/github/KK1182112KK/ponomarev-2016-reproduction/blob/master/python/notebook.ipynb)

```bash
pip install numpy scipy matplotlib pytest
python python/run_reporting_audit.py
python -m pytest python/tests -q
```

The runner executes nine declared cases and writes every trajectory CSV, stdout log, parameters, and source/environment metadata to `results/reporting-audit/`. Full raw runs and trajectories are also in the [CI evidence artifact](https://github.com/KK1182112KK/ponomarev-2016-reproduction/actions/runs/34175280499); that artifact expires after 30 days, while committed summaries and the reproducible runner remain in Git.

MATLAB source is retained. Its entry point is `cd matlab; run_all`, but **MATLAB was not executed in this reporting audit**. The audit notebook executes the same reporting runner; the CI test statement refers to the specific tested source snapshot below.

## Verification and provenance

**32 existing Python tests passed** on the [recorded GitHub Actions run](https://github.com/KK1182112KK/ponomarev-2016-reproduction/actions/runs/34175280499), and all nine declared numerical cases completed both locally and in that run. [Raw test log](results/reporting-audit/tests.log) · [Test record](results/reporting-audit/test-results.json) · [Source/environment manifest](results/reporting-audit/provenance.json).

The solver baseline is `a491223`; the audit runner head is `ece0a41`; CI checked out test merge `46affcfa6c402ff571751af001c79ab6bbbd6365`. The solvers and all pre-existing test assertions were unchanged. The report distinguishes numerical refinement, algebraic identities, and replay of an existing control sequence from stronger independent validation.

[Exact numerical contract](docs/METHODS.md) · [Source specification](docs/SPEC.md) · [Reporting standard](docs/REPORTING_STANDARD.md) · [Related reproduction projects](https://github.com/KK1182112KK/krstic-2016-reproduction)

Implementation code remains under the existing [MIT license](LICENSE). Original theory and figures remain attributed to their authors. No paper PDF or author-generated plot is redistributed.
