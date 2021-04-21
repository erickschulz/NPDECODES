- Whole Problem: Links refer to problem sheet instead of to the lecture document
- Problem Description: The function $\phi$ is implemented in the file `exponentialintegrator.cc`
- (7-5.a) Solution: Typo: Therefor -> Therefore
- (7-5.a) Solution: Typo: whence -> hence
- (7-5.c) Solution: Some $y$ are not bold
- (7-5.d): Wrong function signature: `VectorXd ExpEulStep(const VectorXd &y0, const Function &&f, const Jacobian &&df, double h)` -> `VectorXd exponentialEulerStep(const VectorXd &y0, Function &&f, Jacobian &&df, double h)`
- (7-5.d) Solution: Update code snippet
- (7-5.e): Change function signature: `void testExpEulerLogODE(void)` -> `void testExpEulerLogODE()`
- (7-5.e): Function is implemented in the file `exponentialintegrator.cc`
- (7-5.e): Can move hint to problem formulation, as the exact solution is already coded in the template code anyways


# Error table for (7-5.e)
N = 2      Error = 0.00168674
N = 4      Error = 0.000452805  Approximated order = 1.89728
N = 8      Error = 0.000117603  Approximated order = 1.94497
N = 16     Error = 2.99917e-05  Approximated order = 1.97129
N = 32     Error = 7.57469e-06  Approximated order = 1.98531
N = 64     Error = 1.90346e-06  Approximated order = 1.99256
N = 128    Error = 4.77101e-07  Approximated order = 1.99626
N = 256    Error = 1.1943e-07   Approximated order = 1.99812
N = 512    Error = 2.98771e-08  Approximated order = 1.99906
N = 1024   Error = 7.4717e-09   Approximated order = 1.99953
N = 2048   Error = 1.86823e-09  Approximated order = 1.99977
N = 4096   Error = 4.67095e-10  Approximated order = 1.99988
N = 8192   Error = 1.16777e-10  Approximated order = 1.99996
N = 16384  Error = 2.91974e-11  Approximated order = 1.99984
N = 32768  Error = 7.29886e-12  Approximated order = 2.0001

