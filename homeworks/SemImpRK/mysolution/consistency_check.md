# Inconsistencies

## Subproblem d:
- The function signature was updated to
```
template <class Function, class Jacobian>
std::vector<Eigen::VectorXd> SolveRosenbrock(Function &&f, Jacobian &&df,
                                             const Eigen::VectorXd &y0,
                                             unsigned int M, double T);
```

## Subproblem e: 
The function signature was updated to
```
double CvgRosenbrock();
```
- The function is implemented in the file semimprk.cc
- The resulting error table is:
              M        maxerr 
            160     0.00181187
            320    0.000439645
            640    0.000108218
           1280    2.68249e-05
           2560    6.66039e-06
           5120    1.64229e-06
          10240    3.90643e-07

- The estimated convergence rate is : 2.02461


# (Broken) Links to Lecture Document:
- Introduction Box
- Subproblem c exercise, solution
- Subproblem e exercise
- Subproblem f exercise, solution