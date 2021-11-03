# Inconsistencies.
## Subproblem a):
- The function signature was updated to 
```
Eigen::Vector2d PredPrey(Eigen::Vector2d y0, double T, unsigned N);
```
## Subproblem b):
- The function signature was updated to 
```
void SimulatePredPrey();
```
- Hint 1: Solution no longer uses the general RKIntegrator class, but the implementation of subproblem a) instead.
- Hint 2: Added as a unit test for the PredPrey function, so the hint can be dropped.
# (Broken) Links to Lecture Document:
- Introduction 
- Subproblem c solution
