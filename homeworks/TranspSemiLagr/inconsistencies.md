### Inconsistencies

LehrFEM++ implementation of the BETL problem "Problem 7.5 Transport Problem" in the old problem pdf. Since the description of the old problem is based on BETL, all code-related descriptions are inconsistent. This includes:

* Declaration of function "semiLagr_step" in the introduction
* Declaration of function "solverot"  in subproblem a)
* Solution of subproblem a)
* Declaration of function "reaction_step" in subproblem d)
* Hint for subproblem d) describes BETL classes.
* Solution for subproblem d)
* Declaration of function "solvetrp" in subproblem e)
* Hint for subproblem e) 
* Solution for subproblem e)

Furthermore all references in the original problem pdf refer the old lecture notes.