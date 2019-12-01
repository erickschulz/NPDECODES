/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 06.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "solve.h"

#include <iostream>

#include <Eigen/Core>

#include <lf/uscalfe/uscalfe.h>

#include "MyLinearFEElementMatrix.h"
#include "MyLinearLoadVector.h"

namespace ElementMatrixComputation {

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solvePoissonBVP() {
  // Convert tPoissonda function f to a LehrFEM++ mesh function object
  lf::uscalfe::MeshFunctionGlobal mf_f{f};

  // Define the solution vector
  Eigen::VectorXd solution = Eigen::VectorXd::Zero(1);

  //====================
  // Your code goes here
  //====================

  return solution;
}
/* SAM_LISTING_END_2 */

Eigen::VectorXd solveNeumannEq() {
  // Define the solution vector
  Eigen::VectorXd solution;

  //====================
  // Your code goes here
  //====================

  return solution;
}

}  // namespace ElementMatrixComputation
