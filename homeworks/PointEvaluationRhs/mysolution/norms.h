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
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/mesh_hierarchy.h>
#include <lf/refinement/refutils.h>
#include <lf/uscalfe/uscalfe.h>
#include <mgl2/mgl.h>

namespace PointEvaluationRhs {

double computeH1seminormLinearFE(const lf::assemble::DofHandler &dofh,
                                 const Eigen::VectorXd &mu);

double computeL2normLinearFE(const lf::assemble::DofHandler &dofh,
                             const Eigen::VectorXd &mu);

class MassLocalMatrixAssembler {
 private:
 public:
  explicit MassLocalMatrixAssembler() {}
  bool isActive(const lf::mesh::Entity &entity) const { return true; }
  Eigen::MatrixXd Eval(const lf::mesh::Entity &entity) const;
};

}  // namespace PointEvaluationRhs
#endif  // define __NORMS_H
