Compared against the example solution of the Jan 2020 exam:

- Solution of part e.): `..., if we know the initial data u_o(x) satisfies ...` (instead of: `satisfy`)  **(data is plural)**
- Solution of part e.): `The scheme (0.2.7) has a spatial stencil with m = 1 ...` (instead of: `... with of m = 1 ...`) **(fixed)**
- Solution of part e.): `The CFL-condition enforces the ...` (Instead of: `The the CFL-conditions enforces the ...`) **(fixed)**
- Problem statement in  part h.): Instead of a Cauchy problem we now study initial-value problems for (0.2.1) in `a` 1-periodic setting, ... **(fixed)**
- Hint-1 of part h.): The 'Fig. 2' label is overlapping with the text, and the x-tick labels are overlapping with the x-ticks in the figure. **(ignored)**
- Solution of part h.): The problem statement says that the initial content of this method is just a copy of slopelimfluxdiff(), but the actual initial code in slopelimfluxdiff.h is not a clean copy of slopelimfluxdiff(): The first two highlighted lines in the PDF-code, 'sigma[0] = ...' and 'sigma[n - 1] = ...' are not within #if SOLUTION, thus they are a change to the raw slopelimfluxdiff(). I provide an example replacement for this inconsistency at the end of this file. I tested replacing this part and it works fine for me.  **(fixed, see below)**
- Part j.): Here the students are required to determine the order of the SSP Runge-Kutta timestepping method from some listed errors (the solution of this problem). In other problems where we ask students to determine the order, we do not list the solution explicitly in the PDF. I think we should just remove this listing here and state that the order of convergence should be determined from the code output, or we could add a subproblem where the students should write code in extendedmuscl_main.cc to determine the order from the output computationally. **(maybe we add another subtask for the study of convergence)**
- Solution Code in part k.) needs to be updated (changed `alfa` to `alpha`) **(fixed)**


**(done)**
In slopelimfluxdiff.h, replace lines 82-89 with:
```
#if SOLUTION
  // Computation of slopes \Blue{$\sigma_j$}, uses \Blue{$\mu_{-1}=\mu_{n-1}$},
  // \Blue{$\mu_n=\mu_0$}, which amounts to constant extension of states
  // beyond domain of influence \Blue{$[a,b]$} of non-constant intial data. Same
  // technique has been applied in \lref{cpp:fluxdiff}
  sigma[0] = slopes(mu[n - 1], mu[0], mu[1]); // @\Label[line]{slfd:1}@
#else
  sigma[0] = slopes(mu[0], mu[0], mu[1]);
#endif
  for (int j = 1; j < n - 1; ++j)
    sigma[j] = slopes(mu[j - 1], mu[j], mu[j + 1]);
#if SOLUTION
  sigma[n - 1] = slopes(mu[n - 2], mu[n - 1], mu[0]); // @\Label[line]{slfd:2}@
#else
  sigma[n - 1] = slopes(mu[n - 2], mu[n - 1], mu[n - 1]);
#endif
```
