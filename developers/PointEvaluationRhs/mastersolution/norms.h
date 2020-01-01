#ifndef __NORMS_H
#define __NORMS_H

/**
 * @ file norms.h
 * @ brief NPDE homework PointEvaluationRhs code
 * @ author Christian Mitsch
 * @ date 22.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>

namespace PointEvaluationRhs{

double computeH1seminormLinearFE(const lf::assemble::DofHandler &dofh,
                                 const Eigen::VectorXd &mu);

double computeL2normLinearFE(const lf::assemble::DofHandler &dofh,
                             const Eigen::VectorXd &mu);

class MassLocalMatrixAssembler{
private:
public:
  explicit MassLocalMatrixAssembler() = default;
  bool isActive(const lf::mesh::Entity &entity) { return true; }
  Eigen::MatrixXd Eval(const lf::mesh::Entity &entity);
};

} // namespace PointEvaluationRhs
#endif // define __NORMS_H
