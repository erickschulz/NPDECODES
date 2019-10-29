#ifndef __POINTEVALUTION_H
#define __POINTEVALUTION_H
/**
 * @ file pointEvaluation.h
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

std::pair<double, double> normsSolutionPointLoadDirichletBVP(
    const lf::assemble::DofHandler &dofh, Eigen::Vector2d source_point,
     Eigen::VectorXd &sol_vec);

Eigen::Vector2d GlobalInverseTria(Eigen::Matrix<double, 2, 3> mycorners,
                                  Eigen::Vector2d x);

Eigen::Vector2d GlobalInverseQuad(Eigen::Matrix<double, 2, 4> mycorners,
                                  Eigen::Vector2d x);

inline double triaArea(const Eigen::Vector2d a, const Eigen::Vector2d b,
                       const Eigen::Vector2d c);

std::pair<double, double> solveQuadraticEquation(double a, double b, double c);

class DeltaLocalVectorAssembler {
 private:
  Eigen::Vector2d x_0;
  bool already_found;

 public:
  explicit DeltaLocalVectorAssembler(Eigen::Vector2d x)
      : x_0(x), already_found(false) {}
  bool isActive(const lf::mesh::Entity &entity) const {
    return (!already_found);
  }
  Eigen::VectorXd Eval(const lf::mesh::Entity &entity);
};

}  // namespace PointEvaluationRhs
#endif  // define __POINTEVALUATION_H
