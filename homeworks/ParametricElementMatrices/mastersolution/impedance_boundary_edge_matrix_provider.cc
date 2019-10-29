
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
  result = (*ImpedanceBoundaryEdgeMatrixProvider::bd_flags_)(edge);
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
  // get local to global mapping to obtain the coefficients of w
  const lf::assemble::DofHandler &dofh{
      ImpedanceBoundaryEdgeMatrixProvider::lin_Lagr_fe_space_->LocGlobMap()};
  auto glob_indices = dofh.GlobalDofIndices(cell);
  Eigen::Vector2d w;
  for (int i = 0; i < 2; i++) {
    w(i) = ImpedanceBoundaryEdgeMatrixProvider::w_coeff_vec_(glob_indices[i]);
  }

  // set up the computed matrices
  Eigen::MatrixXd m_1(2, 2);
  m_1 << 24, 6, 6, 4;
  m_1 *= w(0) * w(0);

  Eigen::MatrixXd m_2(2, 2);
  m_2 << 6, 4, 4, 6;
  m_2 *= w(0) * w(1);

  Eigen::MatrixXd m_3(2, 2);
  m_3 << 4, 6, 6, 24;
  m_3 *= w(1) * w(1);

  result = m_1 + m_2 + m_3;

  // multiply by length and divide by 120 as derived
  auto geom = cell.Geometry();
  double length = lf::geometry::Volume(*geom);
  result *= length / 120.;
  /* SOLUTION_END */
  return result;
}
/* SAM_LISTING_END_2 */

}  // namespace ParametricElementMatrices
