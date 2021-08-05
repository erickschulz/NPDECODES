# LehrFEM++ demo codes for Chapter 10.2.2 of the Course "Numerical Methods for PDEs"

The codes in this folder correspond to some examples in Section 10.2.2 of the "Numerical Methods for Partial Differential Equations" lecture notes. 

## Internal layer

This example found in the file `layer_main.cc` performs the experiments 10.2.2.15 & 10.2.2.30 from the lecture nodes. The example examines the Convection-Diffusion problem in the limiting case epsilon->0. A manufactured solution  with an internal layer is approximated by 

- A standard Galerkin Solution
- A Upwind Quadrature Solution
- A Streamline-Diffuion (SUPG) Solution

The example examines the three approximations for spurious oscillations along a diagonal crossing the internal layer.

## Convergence
This example found in the file `convergence_main.cc` performs the experiment 10.2.2.31 from the lecture nodes. The example examines the Convection-Diffusion problem in the case epsilon=1. A smooth manufactured solution is approximated by 

- A standard Galerkin Solution
- A Upwind Quadrature Solution
- A Streamline-Diffuion (SUPG) Solution

The example examines the asymptotic convergence of the three approximations in the L2 norm.