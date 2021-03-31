# Changes
- Changed all function names to fit Coding Guidelines
- Replaced Eigen::Vectors with double, since tests etc. only used length 1 vectors.
- Plotting script plot_tangent.py for the last subproblem
- Plotting script creates two additional plots, that make it more visible, how the timestep is adapted in the solution.
- bugfixes

# Inconsistencies.
## Subproblem b):
- The function is no longer EIGEN-based
- Function signature was updated to
```
template <class DiscEvlOp>
double PsiTilde(const DiscEvlOp& Psi, unsigned int p, double h, double y0);
```
## Subproblem c):
- The function signature was updated to 
```
template <class DiscEvlOp>
std::vector<double> OdeIntEqui(const DiscEvlOp& Psi, double T, double y0,int M);
```
- Hint: ```std::vector<Vector>::push_back()``` -> ```std::vector<double>::push_back()```

## Subproblem d):
- The function signature was updated to 
```
double TestCvpExtrapolatedEuler();
```

## Subproblem e):
- The function signature was updated to
```
template <class DiscEvlOp>
std::pair<std::vector<double>, std::vector<double>> OdeIntSsCtrl(
    const DiscEvlOp& Psi, unsigned int p, double y0, double T, double h0,
    double reltol, double abstol, double hmin);
```

## Subproblem f):
- The function Signature changed to
```
std::pair<std::vector<double>, std::vector<double>> SolveTangentIVP();
```
- Plot was updated.
- 2 New plots to explicitly show adaptive timestep.

# (Broken) Links to Lecture Document:
- Introduction Box
- Subproblem a (exercise, solution)
- Subproblem b exercise
- Subproblem c solution
- Subproblem e exercise, Hint
- Subproblem f solution