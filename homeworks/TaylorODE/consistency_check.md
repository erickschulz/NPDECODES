# Inconsistencies.

## Subproblem c):
- The function signature was updated to 
```
std::vector<Eigen::Vector2d> SolvePredPreyTaylor(const PredPreyModel& model,
                                                 double T,
                                                 const Eigen::Vector2d& y0,
                                                 unsigned int M);
```
- the file name should be updated to "taylorode.h"
- Update hint: They should not update the function signature, otherwise the testing is broken. Instead 
mention that they should implement the missing functions for f,df,d2f in the helper class PredPreyModel.
- Code implemented in new function -> reference to code probably broken.

## Subproblem d):
- The function signature was updated to 
```
double TestCvgTaylorMethod();
```
- Code implemented in new function -> reference to code probably broken.

# Links to old Lecture Document:
- Introduction Box
- Subproblem a hint
- Subproblem b hint2, hint3, 
- Subproblem c exercise

