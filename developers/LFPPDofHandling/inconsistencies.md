Problem 2-9: Handling DoFs

- Link to GITLAB folder in introduction box is not working (new folder name)->
"homeworks/LFPPDofHandling" did not exist on "master"->Folder is now called `HandlingDOFs` **(fixed)**
- All other links to GITLAB also require update (due to change of filenames) **(fixed)**
- Time hints missing for all parts except 2-9.b) **(fixed)**
- function signatures of `NoLocalDofs` and `NoInteriorDofs` need to be updated (now `NumLocalDofs`and `NumInteriorDofs`) **(fixed)**
- typo after the solution of 2-9.b): Extra "The" at the beginning **(fixed)**
- Filename in the beginning of 2-9.b) wrong: Should be `handlingdofs_main.cc` **(fixed)**
- in 2-9.c), maybe add that the function `countEntityDofs`should be adapted in `handlingdofs.cc` **(fixed)**
- in the solution of 2-9.c) two possible solution approaches are listed, maybe add that the following code is approach 2. **(fixed)**
- solution code in 2-9.c) needs to be updated (PDF) **(fixed)**
- In the header of the solution code in 2-9.c), the reference to "... Sub-problem (2-9.c):" is missing (not consistent with other headers in this problem) **(fixed)**
- solution code in 2-9.d) needs to be updated (PDF) **(fixed)**
- In the solution of 2-9.d) it is stated that we loop from `0` to `Ndofs - 1`, but actually we have a for-each loop. Clearer would be something like: "In a loop over all entities of codimensions 1 and 2 in the given mesh, one checks the boundary flag for each entity. If the flag is true, we increment by the number of DOFs associated with this boundary entity." **(fixed)**
- typo in Hint-1 for 2-9.e): "The integral of `bayrcentric` coordinate ..." **(fixed)**
- solution of 2-9.e): Missing brackets after `DofHandler::GlobalDofIndices` -> `DofHandler::GlobalDofIndices()` **(fixed)**
- solution code in 2-9.e) needs to be updated (PDF) **(fixed)**
- solution code in 2-9.f) needs to be updated (PDF) **(fixed)**
- In the header of the solution code in 2-9.f), a ":" is missing after "... Sub-problem (2-9.f)`:`" **(fixed)**
- In the solution of 2-9.g), missing comma after "Nevertheless" **(fixed)**
- solution code in 2-9.g) needs to be updated (PDF) **(fixed)**
