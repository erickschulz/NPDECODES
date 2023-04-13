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
  std::string meshfile = CURRENT_SOURCE_DIR "/../meshes/square.msh";
  if (argc > 1) meshfile = argv[1];
  std::string vtkfile = "graphminsurf.vtk";
  MinimalGraphSurface::graphMinSurfVis(meshfile,vtkfile);
  return 0;
}
