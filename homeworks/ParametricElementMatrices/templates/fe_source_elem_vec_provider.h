/**
 * @file fe_source_elem_vec_provider.cc
 * @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans
 * @date 14/03/2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>

#include <iostream>
#include <stdexcept>

namespace ParametricElementMatrices {

class FESourceElemVecProvider {
 public:
  using elem_vec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  using ElemVec = const elem_vec_t;

  bool isActive(const lf::mesh::Entity &cell) { return true; }

  FESourceElemVecProvider(
      std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>>
          lin_Lagr_fe_space,
      Eigen::VectorXd w_coeff_vec);

  ElemVec Eval(const lf::mesh::Entity &cell);

 private:
  std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> lin_Lagr_fe_space_;
  Eigen::VectorXd w_coeff_vec_;
};

}  // namespace ParametricElementMatrices
