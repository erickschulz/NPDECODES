/**
 * @file local_laplace_qfe.h
 * @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans
 * @date 27/03/2019
 * @copyright Developed at ETH Zurich
 */

#ifndef LOC_LAP_QFE
#define LOC_LAP_QFE

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>

namespace DebuggingFEM {

class LocalLaplaceQFE1 {
 public:
  bool isActive(const lf::mesh::Entity &cell) { return true; }
  Eigen::Matrix<double, 6, 6> Eval(const lf::mesh::Entity &cell);
};

class LocalLaplaceQFE2 {
 public:
  bool isActive(const lf::mesh::Entity &cell) { return true; }
  Eigen::Matrix<double, 6, 6> Eval(const lf::mesh::Entity &cell);
};

class LocalLaplaceQFE3 {
 public:
  bool isActive(const lf::mesh::Entity &cell) { return true; }
  Eigen::Matrix<double, 6, 6> Eval(const lf::mesh::Entity &cell);
};

}  // namespace DebuggingFEM

#endif
