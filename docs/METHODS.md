# Numerical methods — actual baseline implementation

This document supersedes earlier planned/unchecked descriptions. Solver source is unchanged from baseline `a491223df3ebc583e967fa81e4625d18afaf5824`. [The report](REPRODUCTION_REPORT.md) maps operations to paper equations and measured results.

A advances source (73), B advances (75), C advances (89), using **Forward Euler**. A/C reconstruct the spatial predictor at each control update; B evaluates (86)–(88). No independently stabilized target generates the physical state. This does not make the discretization exact.

All three `_lookup_u` helpers use nearest-index rounding, `round(t_query/dt)+N_h`, followed by array-bound clamping. `N_h=round(h/dt)`. This is **not** the timestamp/arrival-time lookup of the separate Krstic project or general zero-order-hold lookup. At h=pi/4 and dt=0.01 the physical grid lag is 79 steps; the predictor spans pi/4 with ds=(pi/4)/79.

A/B use trapezoidal history integrals. Before the new command is stored, the upper endpoint refers to a preallocated zero slot; A recomputes its physical integral after storing the command. A discrete quadrature endpoint has nonzero weight even though one point has zero measure in the continuous integral. This convention is disclosed, not silently repaired.

A advances xi and beta by Euler; f' uses central differences at epsilon=1e-7. If abs(B)<1e-10, code returns zero control. This guard is preserved and its activation count is not logged. C advances xi and beta with the current xi Jacobian (91) and gain equations (92)–(93), then applies printed (104). Final controls are held copies; no gain-factor correction is automatic.

Existing tests and assertions are unchanged. The C Euler/RK4 check replays one command sequence through a second plant stepper; it is not an independently closed-loop controller. B cascade equivalence is an algebraic check. Refinement to a finer Euler run is not a certified exact-solution error. MATLAB, pixel fitting to Fig. 1, and a general robustness/RoA proof were not executed in this audit.

Future interpolation, startup-quadrature or guard corrections should be separate changes with before/after evidence. Do not retroactively call this unchanged baseline a fully repaired numerical implementation.
