
# Inconsistencies.
## Subproblem a):
- The notation f: IxD \subset R^N -> R^N is confusing.

## Subrpoblem c):
- Plot is missing.

## Subproblem f):
- The function signature was updated to 
```
Eigen::Vector2d SdirkStep(const Eigen::Vector2d &z0, double h, double gamma);
```
- The current version of the code reuses the LU-Decomposition of A.


## Subproblem g:
- The function signature was updated to 
```
double CvgSDIRK();
```
- The current solution only tabulates the error at final time T and not the maximal error over all timesteps.


# (Broken) Links to Lecture Document:
- Introduction Box
- Subproblem a exercise, hint
- Subproblem b solution
- Subproblem d solution
