
/**
 * @file impedance_boundary_edge_matrix_provider.cc
 * @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans
 * @date 14/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "impedance_boundary_edge_matrix_provider.h"

namespace ParametricElementMatrices {

/**
 * @brief initialize ImpedanceBoundaryEdgeMatrixProvider for \int_{\boundary
 * \Omega} w(x)^2 u(x) v(x) dx
 * @param lin_Lagr_fe_space discretization of computational domain
 *        w_coeff_vec coefficients of w interpolated on lin_Lagr_fe_space
 */
/* SAM_LISTING_BEGIN_0 */
ImpedanceBoundaryEdgeMatrixProvider::ImpedanceBoundaryEdgeMatrixProvider(
    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>>
        lin_Lagr_fe_space,
    Eigen::VectorXd w_coeff_vec) {
  ImpedanceBoundaryEdgeMatrixProvider::lin_Lagr_fe_space_ = lin_Lagr_fe_space;
  ImpedanceBoundaryEdgeMatrixProvider::w_coeff_vec_ = w_coeff_vec;
  // initialize boundary flags for isActive() function
  auto mesh = lin_Lagr_fe_space->Mesh();
  ImpedanceBoundaryEdgeMatrixProvider::bd_flags_ =
      std::make_shared<lf::mesh::utils::CodimMeshDataSet<bool>>(
          lf::mesh::utils::flagEntitiesOnBoundary(mesh, 1));
}
/* SAM_LISTING_END_0 */

/**
 * @brief true if edge is on the boundary
 * @param edge edge to be tested for
 */
/* SAM_LISTING_BEGIN_1 */
bool ImpedanceBoundaryEdgeMatrixProvider::isActive(
    const lf::mesh::Entity &edge) {
  bool result;
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return result;
}
/* SAM_LISTING_END_1 */

/**
 * @brief compute element matrix for \int_{\boundary \Omega} w(x)^2 u(x) v(x) dx
 * @param cell edge on the boundary
 */
/* SAM_LISTING_BEGIN_2 */
ImpedanceBoundaryEdgeMatrixProvider::ElemMat
ImpedanceBoundaryEdgeMatrixProvider::Eval(const lf::mesh::Entity &cell) {
  Eigen::MatrixXd result(2, 2);
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return result;
}
/* SAM_LISTING_END_2 */

}  // namespace ParametricElementMatrices
