# Specification: Ponomarev (2016)

## Reference
- **Authors**: Anton Ponomarev
- **Title**: Nonlinear Predictor Feedback for Input-Affine Systems with Distributed Input Delays
- **Venue**: IEEE Transactions on Automatic Control
- **Year**: 2016
- **DOI**: arXiv:1601.00098v1 (accepted to IEEE TAC, Oct 2015)

## Problem Statement

This paper extends predictor feedback methodology to control-affine nonlinear systems with distributed input delays. The key contribution is a state transformation that maps the delayed system to a delay-free form, enabling standard stabilization techniques. Three examples demonstrate the approach: scalar analytical, explicit prediction, and numerical prediction for an inverted pendulum.

## General System Model (Eq. 4)

```
dx/dt = f(x) + B0(x)*u(t) + B1(x)*u(t-h) + ∫_{-h}^{0} Bint(θ,x)*u(t+θ) dθ
```

where x ∈ R^n, u ∈ R^m, h > 0.

### Predictor Transformation (Eq. 31-34)

y(t) = Y(x(t), ut) where Y(x, ϕ) = ξ(h) with:
```
ξ'(s) = f(ξ(s)) + B1(ξ(s))*ϕ(s-h) + ∫_{-h}^{-s} Bint(θ,ξ(s))*ϕ(s+θ) dθ
ξ(0) = x
```
Solve for s ∈ [0, h].

### Transformed System (Theorem 1, Eq. 39-45)

```
dy/dt = f(y) + B(y, ut)*u(t)
```

where B(y, ϕ) = B1(ξ(h)) + β(h) with:
```
β'(s) = Ã(s, ξ(s), ϕ)*β(s) + Bint(-s, ξ(s))
β(0) = B0(ξ(0))
ξ'(s) = [predictor ODE, integrated from ξ(h)=y backwards]
```

and
```
Ã(s, ξ, ϕ) = A(ξ) + Σ_{i=1}^{m} Bi_1(ξ)*ϕi(s-h) + Σ_{i=1}^{m} ∫_{-h}^{-s} Bi_int(θ,ξ)*ϕi(s+θ) dθ
```

### Feedback Law (Corollary 1, Eq. 71)

```
κ(y, ϕ) = -k * B^T(y, ϕ) * ∇v0(y)
```

---

## Example A: Scalar Case (Section VI-A, Eq. 73-74)

### System
```
dx/dt = f(x) + b0*u(t) + b1*u(t-h) + ∫_{-h}^{0} bint(θ)*u(t+θ) dθ
```
Scalar system (n=1, m=1) with same-sign coefficients: b0 > 0, b1 > 0, bint > 0.

### Predictor Transformed System
```
dy/dt = f(y) + B(y,ϕ)*u(t)
```
where B(y,ϕ) > 0 (positive scalar, bounded away from zero).

### Feedback (Eq. 74)
```
κ(y, ϕ) = -(f(y) + y) / B(y, ϕ)
```
This yields dy/dt = -y (exponential decay).

### Test Parameters (not specified in paper — chosen for demonstration)
| Parameter | Value | Description |
|-----------|-------|-------------|
| f(x) | sin(x) | Nonlinear open-loop dynamics |
| b0 | 1.0 | Undelayed input coefficient |
| b1 | 0.5 | Delayed input coefficient |
| bint(θ) | 1.0 | Distributed delay kernel (constant) |
| h | 0.5 s | Input delay |
| x(0) | 1.0 | Initial state |
| u(θ) | 0 for θ ∈ [-h, 0] | Input history |
| t_end | 10 s | Simulation time |
| dt | 0.001 s | Euler step size |

---

## Example B: Explicit Prediction (Section VI-B, Eq. 75-88)

### System (Eq. 75)
```
dx1/dt = x2^2 + u(t-h)
dx2/dt = x2 + u(t)
```

### Predictor Transformation (Eq. 76-77)
```
y1 = x1 + (e^(2h) - 1)/2 * x2^2 + ∫_{-h}^{0} u(t+θ) dθ
y2 = e^h * x2
```

### Transformed System (Eq. 78)
```
dy1/dt = y2^2 + (1 + (e^h - e^(-h))*y2) * u(t)
dy2/dt = y2 + e^h * u(t)
```

### Cascade Transformation (Eq. 79-81)
```
z1 = y1 - e^(-h)*y2 + (e^(-2h) - 1)/2 * y2^2
z2 = e^(-h) * y2
u = -2*z2 + ũ
```

### Cascade System (Eq. 82)
```
dz1/dt = z2^2 - z2
dz2/dt = -z2 + ũ
```

### Lyapunov Function (Eq. 83)
```
V(z) = (z1 + z2*(z2-2)/2)^2 + z2^2
```

### Feedback (Eq. 84-85)
```
dV/dt = -2*z2^2 + ∂V/∂z2 * ũ
ũ = -∂V/∂z2
```

### Overall Simplified Transformation (Eq. 86-88)
```
z1 = x1 - x2 + ∫_{-h}^{0} u(t+θ) dθ
z2 = x2
u = -2*x2 + ũ
```
where ũ = -∂V/∂z2 with V from Eq. 83.

### Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| h | 1.0 s | Input delay |
| x(0) | [1, 1]^T | Initial state |
| u(θ) | 0 for θ ∈ [-h, 0] | Input history |
| t_end | 15 s | Simulation time |
| dt | 0.001 s | Euler step size |

---

## Example C: Numerical Prediction / Inverted Pendulum (Section VI-C, Eq. 89-104)

### System (Eq. 89)
```
dx1/dt = x2
dx2/dt = sin(x1) + u(t) + u(t-h)
h = π/4
```

### System Matrices (Eq. 90-91)
```
f(x) = [x2; sin(x1)]
B0 = B1 = [0; 1]
Bint = 0 (no distributed delay integral term)

A(x) = [0, 1; cos(x1), 0]
Mf = 1
R = ∞ (global)
```

### β ODE (Eq. 92-93)
Since Bint = 0, the Ã matrix simplifies:
```
β'(s) = A(ξ(s)) * β(s)
β(0) = B0 = [0; 1]
```

### Predictor ODE (Eq. 46-47, specialized)
Since Bint = 0:
```
ξ'(s) = f(ξ(s)) + B1 * u(t + s - h)
       = [ξ2(s); sin(ξ1(s)) + u(t + s - h)]
ξ(0) = x(t)
```

### B Computation (Eq. 40)
```
B(y, ϕ) = B1 + β(h) = [0; 1] + β(h)
```

### Feedback (Eq. 104, via Corollary 1)
Paper states V = I, k = 1 and writes the feedback as:
```
κ(y, ϕ) = -B^T(y, ϕ) * y
```
Note: The general formula Eq. 71 gives κ = -k*B^T*∇v0 = -2k*B^T*V*y.
With V=I, k=1 this would be -2B^T*y. The paper's Eq. 104 omits the
factor of 2, effectively using k=1/2 or absorbing it into the convention.
We implement Eq. 104 as written to match Fig. 1.

### Control Algorithm (5 steps, from paper)
At each time t:
1. Given x(t) and ut, solve predictor ODE for ξ(s), s ∈ [0,h], with ξ(0)=x(t)
2. Set y(t) = ξ(h)
3. Solve β ODE for β(s), s ∈ [0,h], with β(0)=[0;1]
   Note: β ODE uses the SAME ξ(s) trajectory from step 1
4. Compute B(y,ut) = [0;1] + β(h)
5. Apply u(t) = -B^T(y,ut) * y(t)

### Parameters (from paper, Fig. 1)
| Parameter | Value | Description |
|-----------|-------|-------------|
| h | π/4 ≈ 0.7854 s | Input delay |
| x1(0) | {π/2, π, 3π/2} | Initial angle (3 cases in Fig. 1) |
| x2(0) | π/2 | Initial angular velocity |
| u(θ) | 1 for θ ∈ [-h, 0] | Input history (constant 1) |
| t_end | 10 s | Simulation time |
| dt | 0.01 s | Euler integration step (per paper) |
| V | I (2×2 identity) | Lyapunov matrix |
| k | 1 | Feedback gain |

---

## Figures to Reproduce
| Figure | Description | Status |
|--------|-------------|--------|
| Fig. 1 (left) | Control input u(t) for Example C, x1(0) ∈ {π/2, π, 3π/2} | [ ] |
| Fig. 1 (right) | State x1(t) for Example C, same 3 initial conditions | [ ] |
| Fig. 2 (new) | Example A: scalar state & control trajectories | [ ] |
| Fig. 3 (new) | Example B: state trajectories x1,x2 and cascade z1,z2 | [ ] |
| Fig. 4 (new) | Example C: compensated vs uncompensated comparison | [ ] |
| Fig. 5 (new) | Example C: parameter sweep over h ∈ [0.1, π/2] | [ ] |

## Success Criteria
1. **Example C matches Fig. 1**: x1 converges from {π/2, π, 3π/2} to 0; u(t) profile matches paper qualitatively (u → 0, initial spike to ≈ -8 or -9)
2. **All states converge**: |x(t_end)| < 1e-2 for all examples
3. **Euler vs ode45 agreement**: max|x_euler - x_ode45| < 0.05 for Example C with same parameters
4. **Example B transformation equivalence**: z-transform (86-88) matches y-transform (76-77) to machine precision
5. **Convergence under dt refinement**: Euler error decreases as O(dt) for Example C
6. **Lyapunov decay**: V(z(t)) and v0(y(t)) decrease monotonically (after initial transient) in Examples B and C
