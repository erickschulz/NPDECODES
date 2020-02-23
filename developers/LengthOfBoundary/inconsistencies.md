Compared against `https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE19_Homework_withSolutions.pdf`:


- in the problem statement of 2-7.a) the function signature `volumeOfDomain` needs to be updated (mesh input is now a pointer & not consistent with solution code), and the reference to `mesh` in the description is wrong (should be `mesh_p`)
- Example Code in the solution of subproblem 2-7.a) needs to be updated (mesh input is now a pointer)
- solution of 2-7.b): Example Code from LF++ (`CountNumSuperEntities`) needs to be updated: still old function name there (`CountNoSuperEntities` vs. `CountNumSuperEntities`). Also the reference to this function in the description is wrong
- typo in solution of 2-7.b): adjcnt -> adjacent
- in the problem statement of 2-7.c) the function signature `lengthOfBoundary`needs to be updated (mesh input is now a pointer & not consistent with solution code), and the reference to `mesh` in the description is wrong (should be `mesh_p`)
- Example Code in the solution of subproblem 2-7.c) needs to be updated
- in the problem statement of 2-7.d) the function signature `measureDomain` needs to be updated (variable name changed), and in the description `msh_file_name` needs to be replaced with `filename`
- in the problem statement of 2-7.d) it is required to edit the `main()` routine, maybe add that this can be found in `boundarylength_main.cc`
- Example Code in the solution of subproblem 2-7.d) needs to be updated: Replace `measureDomain` code and (maybe) add the solution found in the file `boundarylength_main.cc` here too (since it is also asked for in the description)