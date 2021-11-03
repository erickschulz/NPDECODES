# Inconsistencies

## Between Subproblem b and c:
- Butcher Matrix is implemented in gradientflow.cc


## Subproblem e: 
- Since Subrpoblem c) and d) seem to be dropped, this subproblem probably needs an introduction for the function ```ComputeStages```. This function has the new signature
```
template <typename Func, typename Jac>
std::array<Eigen::VectorXd, 5> ComputeStages(Func &&f, Jac &&df,
                                             const Eigen::VectorXd &y, double h,
                                             double rtol = 1E-6,
                                             double atol = 1E-8);
```
- The function signature was updated to
```
template <typename Func, typename Jac>
Eigen::VectorXd DiscEvolSDIRK(Func &&f, Jac &&df, const Eigen::VectorXd &y,
                              double h, double rtol = 1E-6,
                              double atol = 1E-8);
```

## Subproblem h:
- The function signature was updated to
```
std::vector<Eigen::VectorXd> SolveGradientFlow(const Eigen::VectorXd &d,
                                               double lambda,
                                               const Eigen::VectorXd &y0,
                                               double T, unsigned int M);
```
-  N->M ;  ```y```->```y0```
- ```computeStages()``` -> ```ComputeStages()```


# (Broken) Links to Lecture Document:
- Introduction Box
- (Subproblem c exercise, solution) 
- (Subproblem d exercise)
- Subproblem e solution
- Subproblem g hint, solution
