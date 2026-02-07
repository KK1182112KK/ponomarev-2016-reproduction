# Numerical Methods & Implementation Notes

## ODE Solver

### Primary: Forward Euler (per paper)

The paper explicitly states "Euler's approximation with time step of 0.01" for Example C.
We use forward Euler as the primary method to match the paper's results.

```
x(t+dt) = x(t) + dt * dx/dt(t)
```

### Secondary: MATLAB ode45 (cross-validation)

For validation, we also solve with ode45 (adaptive RK45):
- RelTol = 1e-8, AbsTol = 1e-10
- MaxStep = dt (to properly resolve control updates)

### Stiffness Analysis

**Example A** (scalar): Non-stiff. Open-loop dynamics f(x)=sin(x) has bounded derivative.
**Example B** (2-state): Mildly stiff. The dx2/dt = x2 + u(t) term has eigenvalue +1 (unstable open-loop), but the feedback stabilizes it. Use small dt for Euler.
**Example C** (pendulum): Non-stiff. Eigenvalues of A(x) = [0,1; cos(x1),0] have magnitude ≤ 1. Euler with dt=0.01 is adequate.

## Discretization

### Time Stepping

| Example | dt (Euler) | Justification |
|---------|-----------|---------------|
| A | 0.001 s | Small for accuracy with distributed delay integral |
| B | 0.001 s | Small for stability (unstable open-loop mode) |
| C | 0.01 s | Per paper specification |

### Control History Management

All three examples require storing the input history u(t+θ) for θ ∈ [-h, 0].

**Storage**: Array of past control values, indexed by time step.
- At each step, append new u(t) and trim history older than t-2h (keep buffer beyond t-h for interpolation)
- For Example A: also need the distributed delay integral ∫_{-h}^{0} bint(θ)*u(t+θ) dθ, computed via trapezoidal rule over stored history

### Predictor ODE Integration

At each time step, the predictor ODE must be solved in "spatial" variable s ∈ [0, h]:

**Example A**: Predictor ODE with distributed delay kernel
```
ξ'(s) = f(ξ(s)) + b1 * u(t+s-h) + ∫_{-h}^{-s} bint(θ) * u(t+s+θ) dθ
ξ(0) = x(t)
```
Integrate with Euler in s, using N_pred = ceil(h/dt) sub-steps.

**Example B**: Analytical predictor — no ODE to solve.

**Example C**: Two coupled ODEs (predictor + β):
```
ξ'(s) = [ξ2; sin(ξ1) + u(t+s-h)],   ξ(0) = x(t)
β'(s) = A(ξ(s)) * β(s),               β(0) = [0;1]
```
Both integrated simultaneously with Euler in s.
N_pred sub-steps in [0, h], with ds = h / N_pred.
Use N_pred = ceil(h / dt) to maintain consistent resolution.

## Numerical Pitfalls

### 1. Input History Interpolation (All Examples)

The predictor ODE needs u(t+s-h) at arbitrary s values. Since we only store u at discrete time steps, interpolation is needed.

**Workaround**: Use piecewise-constant (zero-order hold) interpolation. Since the paper uses Euler with fixed dt, the control signal is piecewise constant anyway.

### 2. Initial Control Discontinuity (Example C)

At t=0, u jumps from the initial history value (u=1) to the feedback-computed value (≈ -8 to -9). This creates a sharp transient in ξ(s) during the first few time steps.

**Workaround**: No special treatment needed — Euler handles this if dt is small enough. The paper reproduces this behavior in Fig. 1.

### 3. Distributed Delay Quadrature (Example A)

Computing ∫_{-h}^{0} bint(θ)*u(t+θ) dθ from discrete history requires numerical quadrature.

**Workaround**: Trapezoidal rule over stored history values. With constant bint and N_h = h/dt history points:
`bint * dt * (u_0/2 + u_1 + ... + u_{N-1} + u_N/2)` (standard trapezoidal weighting).

### 4. β ODE Coupled to ξ Trajectory (Example C)

The β ODE (92) depends on ξ(s) at each s-step. Both must be integrated together.

**Workaround**: Integrate ξ and β simultaneously in the same Euler loop over s ∈ [0,h].

### 5. Example B: Cascade Feedback Complexity

The feedback for Example B involves ∂V/∂z2 which expands to:
```
∂V/∂z2 = 2*(z1 + z2*(z2-2)/2)*(z2-1) + 2*z2
```
This is a moderately complex expression. Sign errors here will cause instability.

**Workaround**: Implement symbolically, verify by finite differences.

### 6. Control History Growing Unbounded

All examples accumulate control history over time.

**Workaround**: Trim history to [t-2h, t] window. Only [t-h, t] is needed for the predictor, but keep a small buffer for interpolation.

## Validation Strategy

### Three-Tier Testing

**Tier 1: Analytical Comparison**
- Example B: Verify the analytical predictor (76-77) produces the same y as the numerical predictor ODE (32-34) for the same system
- Example C at t=0: With u(θ)=1 for θ∈[-h,0], the predictor ODE can be solved semi-analytically as a forced pendulum
- Example A: With f(x)=0, the predictor is linear and has a closed-form solution

**Tier 2: Convergence Order**
- Example C: Run with dt ∈ {0.04, 0.02, 0.01, 0.005, 0.0025}
- Error metric: max|x_dt - x_ref| where x_ref is ode45 with tight tolerances
- Expected order: O(dt) for forward Euler

**Tier 3: Cross-Method Validation**
- Example C: Euler vs ode45 agreement at dt=0.01
- Example B: Explicit predictor (76-77) vs numerical predictor ODE (32-34)
- Example B: Simplified transformation (86-88) vs full transformation (76-77)→(79-81)

### Lyapunov Verification
- Example B: Monitor V(z(t)) from Eq. 83 — must decrease after transient
- Example C: Monitor v0(y) = ||y||^2 — must decrease after transient
- Example C: Monitor Lyapunov-Krasovskii functional (65) — must decrease

## Performance Notes

- Example A: ~1 second for 10s simulation (scalar, simple)
- Example B: ~1 second (analytical predictor, no ODE in predictor)
- Example C: ~5-10 seconds for 10s simulation
  - At each of 1000 time steps, solve 2 ODEs over ~79 sub-steps
  - Total: ~79,000 Euler steps for predictor + β
- Memory: Control history requires ~1000 * m doubles per example (negligible)
- Figure generation: ~30 seconds for full suite (includes parameter sweep)
