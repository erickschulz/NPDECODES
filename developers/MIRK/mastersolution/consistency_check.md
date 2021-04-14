# Inconsistencies.

## Subrpoblem e):
- The Function signature was updated to
```
template <class Func, class Jac>
Eigen::VectorXd Newton2Steps(Func &&F, Jac &&DF, Eigen::VectorXd z);
```
Remark: The function signature in the current version (with a ```void``` return) resulted in the following problem: If an ```Eigen::Vector3d``` is passed to the function,  the compiler seems not to be able to convert it to an ```Eigen::VectorXd&```. This version of the funciton signature avoids this problem.


## Subproblem f):
- The function newton2steps is now called Newton2Steps (2x)

- The function signature was updated to 
```
template <class Func, class Jac>
double MIRKStep(Func &&f, Jac &&df, double y0, double h);
```


## Subproblem g:
- The function signature was updated to 
```
template <class Func, class Jac>
double MIRKSolve(Func &&f, Jac &&df, double y0, double T, unsigned int M);
```


## Subproblem h:
- n->M (for y_M and M=4,...,512)
- The last line of the error table was missing:
    - 512, 1.55741, 2.5743e-06


# (Broken) Links to Lecture Document:
- Subproblem a: exercise, hint
- Subproblem b: exercise, hint,solution
- Subproblem d: solution 