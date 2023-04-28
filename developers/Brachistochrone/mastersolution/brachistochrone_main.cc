/**
 * @ file
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Core>
#include <iostream>

#include "brachistochrone.h"

int main(int /*argc*/, char** /*argv*/) {
  // Track progress of iteration
  Brachistochrone::tabiterr();
  // Convergence study
  Brachistochrone::brachistochrone_cvg("brachistochrone.m");
  Brachistochrone::iteration_test("iterates.m");
  return 0;
}
