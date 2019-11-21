/**
 * @ file norms.cc
 * @ brief NPDE homework PointEvaluationRhs code
 * @ author Christian Mitsch
 * @ date 22.03.2019
 * @ copyright Developed at ETH Zurich
 */
#include "norms.h"
#include <cmath>

namespace PointEvaluationRhs
{

/* SAM_LISTING_BEGIN_1 */
double computeL2normLinearFE(const lf::assemble::DofHandler &dofh,
                             const Eigen::VectorXd &mu)
{
  double result = 0.0;
  //====================
  // Your code goes here
  //====================
  return result;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double computeH1seminormLinearFE(const lf::assemble::DofHandler &dofh,
                                 const Eigen::VectorXd &mu)
{
  // calculate stiffness matrix by using the already existing local assembler
  // LinearFELaplaceElementMatrix
  double result = 0.0;
  //====================
  // Your code goes here
  //====================
  return result;
}
/* SAM_LISTING_END_2 */

Eigen::MatrixXd MassLocalMatrixAssembler::Eval(
    const lf::mesh::Entity &entity)
{
  Eigen::MatrixXd result;
  //====================
  // Your code goes here
  //====================
  return result;
}

} // namespace PointEvaluationRhs
