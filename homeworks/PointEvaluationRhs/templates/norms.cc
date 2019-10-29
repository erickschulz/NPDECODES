/**
 * @ file norms.cc
 * @ brief NPDE homework PointEvaluationRhs code
 * @ author Christian Mitsch
 * @ date 22.03.2019
 * @ copyright Developed at ETH Zurich
 */
#include "norms.h"

namespace PointEvaluationRhs {

double computeL2normLinearFE(const lf::assemble::DofHandler &dofh,
                             const Eigen::VectorXd &mu) {

  double result = 0.0;
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return result;
}

double computeH1seminormLinearFE(const lf::assemble::DofHandler &dofh,
                                 const Eigen::VectorXd &mu) {
  // calculate stiffness matrix by using the already existing local assembler
  // LinearFELaplaceElementMatrix
  double result = 0.0;
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return result;
}

Eigen::MatrixXd MassLocalMatrixAssembler::Eval(const lf::mesh::Entity &entity) const{
  Eigen::MatrixXd result;
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return result;
}

}  // namespace PointEvaluationRhs
