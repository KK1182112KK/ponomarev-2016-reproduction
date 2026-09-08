# Source specification — Ponomarev 2016

This specification supersedes the initial draft's unchecked success criteria. Equation numbers refer to **arXiv:1601.00098v1**; see the [source-based report](REPRODUCTION_REPORT.md) for the equations and complete mappings.

| Case | Physical equations | Feedback / predictor | Source versus repository conditions |
|---|---|---|---|
| A | (73) | (74), (33)–(34), (40)–(42) | The paper specifies a scalar class with same-sign coefficients. Sine nonlinearity, b0=1, b1=0.5, bint=1, h=0.5, x0=1, zero history, T=10 and dt are repository choices. |
| B | (75) | (83), (85)–(88) | The paper supplies the analytical transform/design; h=1, x0=[1,1], zero history, T=15 and dt are repository choices. No separate published B figure. |
| C | (89) | (46)–(47), (40), (92)–(93), printed (104) | Fig. 1 supplies the three angles, x2(0)=pi/2 and input history 1. h=pi/4 is in (89), Euler 0.01 in text after (104). Additional steps are audit experiments. |

Source-based algebra and numerical implementation are separate. Equation (71) with (98) has a normalization difference from printed (104) for the stated V,k; code retains (104). A finite terminal norm is not Theorem 2's bound (64). Actual pendulum terminal tests use 0.1; the stronger old README claim of 0.01 for every case is not supported.

[All measured cases](../results/reporting-audit/summary.csv) · [Implementation conventions](METHODS.md) · [Reporting standard](REPORTING_STANDARD.md)
