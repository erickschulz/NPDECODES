/**
 * @ file minimalgraphsurface_main.cc
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author R. Hiptmair & W. Tonnon
 * @ date April 2023
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Core>
#include <iostream>

#include "minimalgraphsurface.h"

int main(int argc, char** argv) {
  std::string meshfile = "/../meshes/square.msh";
  if (argc > 1) meshfile = argv[1];
  std::string vtkfile = "graphminsurf.vtk";
  MinimalGraphSurface::graphMinSurfVis(meshfile,vtkfile);
  return 0;
}

/*
$ git push --set-upstream origin MinimalGraphSurface
Enumerating objects: 13, done.
Counting objects: 100% (13/13), done.
Delta compression using up to 8 threads
Compressing objects: 100% (10/10), done.
Writing objects: 100% (11/11), 1.40 KiB | 287.00 KiB/s, done.
Total 11 (delta 3), reused 1 (delta 0), pack-reused 0
remote: Resolving deltas: 100% (3/3), completed with 2 local objects.
remote:
remote: Create a pull request for 'MinimalGraphSurface' on GitHub by visiting:
remote: https://github.com/erickschulz/NPDECODES/pull/new/MinimalGraphSurface
remote:
To https://github.com/erickschulz/NPDECODES.git
 * [new branch]      MinimalGraphSurface -> MinimalGraphSurface
branch 'MinimalGraphSurface' set up to track 'origin/MinimalGraphSurface'.
*/
