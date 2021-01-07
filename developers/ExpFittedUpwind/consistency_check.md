### Inconsistency Check
* LehrFEM++ implementation of the BETL problem "Problem 7.4 Exponentially fitted upwind scheme" in the old problem pdf. Since the description of the old problem is based on BETL, all code-related descriptions are inconsistent. These are the exercise descriptions and solutions for the subproblems j-m. And the solution to subproblem n.
* All references in the original problem pdf refer the old lecture notes.
* Hint 2 for subproblem h): The hint should be changed from "For Psi=0 you should get the element matrix A_K [...]" to "For Psi=0 you should get the element matrix -A_K [...]" (Reason: In this case j = -grad u)
* Solution to subproblem h): The solution formula at the end should contain 
  a matrix multiplication with the diagonal matrix diag(exp(-Psi(a^1_K)),exp(-Psi(a^2_K))exp(-Psi(a^2_K))), instead of diag(exp(-Psi(a^1_K)),exp(-Psi(a^1_K))exp(-Psi(a^1_K))).
* Solution to subproblem n): Code contains an additional python scripts which plots the L2 error. The
corresponding plot could be added to the solution.