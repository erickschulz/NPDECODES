
/**
 * @file impedance_boundary_edge_matrix_provider.h
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

class ImpedanceBoundaryEdgeMatrixProvider {
 public:
  using elem_mat_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  using ElemMat = const elem_mat_t;
  ImpedanceBoundaryEdgeMatrixProvider(
      std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>>
          lin_Lagr_fe_space,
      Eigen::VectorXd w_coeff_vec);
  bool isActive(const lf::mesh::Entity &edge);
  ElemMat Eval(const lf::mesh::Entity &cell);

 private:
  std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> lin_Lagr_fe_space_;
  Eigen::VectorXd w_coeff_vec_;
  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<bool>> bd_flags_;
};

}  // namespace ParametricElementMatrices
