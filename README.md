# Nonlinear Predictor Feedback for Input-Affine Systems with Distributed Input Delays
### Reproduction of Ponomarev (IEEE TAC, 2016)

[![MATLAB Tests](https://github.com/YOUR_USERNAME/ponomarev-2016/actions/workflows/matlab-ci.yml/badge.svg)](https://github.com/YOUR_USERNAME/ponomarev-2016/actions/workflows/matlab-ci.yml)
[![Python Tests](https://github.com/YOUR_USERNAME/ponomarev-2016/actions/workflows/python-ci.yml/badge.svg)](https://github.com/YOUR_USERNAME/ponomarev-2016/actions/workflows/python-ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Overview

This paper extends predictor feedback methodology to control-affine nonlinear systems with distributed input delays. A state transformation maps the delayed system to a delay-free form, enabling standard stabilization techniques. Three numerical examples demonstrate the approach: a scalar system with analytical predictor, a 2-state system with explicit cascade transformation, and an inverted pendulum with numerical predictor integration.

## Key Results

<!-- Hero figures will be added after MATLAB generates them -->
| Example A (Scalar) | Example C (Inverted Pendulum) |
|:---:|:---:|
| ![Example A](matlab/results/fig2_example_a.png) | ![Example C](matlab/results/fig1_paper_figure1.png) |

**Figure 1 (right)** reproduces the paper's Figure 1: predictor feedback stabilizes the inverted pendulum from initial angles x1(0) in {pi/2, pi, 3pi/2} with delay h = pi/4.

## Methods

The predictor transformation Y(x, ut) maps the delayed input system to a delay-free form:

- **Predictor ODE**: Solve xi'(s) = f(xi) + B1*phi(s-h) + integral terms, s in [0, h]
- **Effective gain**: B(y, phi) = B1 + beta(h), where beta solves a coupled ODE
- **Feedback**: kappa = -k * B^T * grad(v0(y))

All three examples use Forward Euler integration (per paper specification, dt = 0.01 for Example C).

## Quick Start

### MATLAB
```matlab
git clone https://github.com/YOUR_USERNAME/ponomarev-2016.git
cd ponomarev-2016/matlab
run_all          % Simulate all 3 examples
run_all('fig')   % Generate figures
run_all('test')  % Run validation suite (18 tests)
run_all('all')   % Everything
```

### Python
```bash
cd ponomarev-2016
pip install -r python/requirements.txt
python python/run_all.py             # Simulate all examples
python python/run_all.py --mode fig  # Generate figures
pytest python/tests/ -v              # Run tests
```

## Validation

| Criterion | Result |
|-----------|--------|
| Example C convergence (3 ICs) | x(t_end) < 1e-2 for all |
| Euler vs RK4 agreement | max error < 0.1 at dt=0.01 |
| Convergence order (Euler) | O(dt) verified across 4 refinements |
| Lyapunov V(z) decay | Non-increasing after transient (>80% steps) |
| Cascade identity (Eq. 86-88 vs 76-77) | Matches to machine precision |
| MATLAB/Python agreement | Cross-validated |

## Examples

### Example A: Scalar Predictor (Eq. 73-74)
- System: dx/dt = sin(x) + u(t) + 0.5*u(t-0.5) + integral(u)
- Feedback: kappa = -(sin(y) + y) / B, yielding dy/dt = -y

### Example B: Explicit Prediction (Eq. 75-88)
- 2-state system with cascade transformation to z-coordinates
- Lyapunov-based feedback with V(z) = (z1 + z2(z2-2)/2)^2 + z2^2

### Example C: Inverted Pendulum (Eq. 89-104)
- System: dx1/dt = x2, dx2/dt = sin(x1) + u(t) + u(t-pi/4)
- Numerical predictor integration at each time step
- Reproduces paper's Figure 1

## Citation

If you use this code, please cite the original paper:

```bibtex
@article{ponomarev2016nonlinear,
  author  = {Anton Ponomarev},
  title   = {Nonlinear Predictor Feedback for Input-Affine Systems with Distributed Input Delays},
  journal = {IEEE Transactions on Automatic Control},
  year    = {2016},
  note    = {arXiv:1601.00098v1}
}
```

## License

This project is licensed under the MIT License -- see [LICENSE](LICENSE) for details.
