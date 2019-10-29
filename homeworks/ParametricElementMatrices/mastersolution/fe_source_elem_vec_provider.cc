/**
 * @file fe_source_elem_vec_provider.cc
 * @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans
 * @date 14/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "fe_source_elem_vec_provider.h"

namespace ParametricElementMatrices {

/**
 * @brief initialize FESourceElemVecProvider for \int_{\Omega} v(x)/(1 + w(x)^2)
 * dx
 * @param lin_Lagr_fe_space discretization of computational domain
 *        w_coeff_vec coefficients of w interpolated on lin_Lagr_fe_space
 */
FESourceElemVecProvider::FESourceElemVecProvider(
    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>>
        lin_Lagr_fe_space,
    Eigen::VectorXd w_coeff_vec) {
  FESourceElemVecProvider::w_coeff_vec_ = w_coeff_vec;
  FESourceElemVecProvider::lin_Lagr_fe_space_ = lin_Lagr_fe_space;
}

/**
 * @brief compute element vector for \int_{\Omega} v(x)/(1 + w(x)^2) dx
 * @param cell edge on the boundary
 */
FESourceElemVecProvider::ElemVec FESourceElemVecProvider::Eval(
    const lf::mesh::Entity &cell) {
  auto geom = cell.Geometry();
  auto ref_element = geom->RefEl();
  Eigen::VectorXd result;
  /* SOLUTION_BEGIN */
  // get local to global map
  const lf::assemble::DofHandler &dofh{
      FESourceElemVecProvider::lin_Lagr_fe_space_->LocGlobMap()};
  auto glob_indices = dofh.GlobalDofIndices(cell);
  // distinguish between quadrilaterals and triangles
  if (ref_element.NumNodes() == 3) {
    /* SAM_LISTING_BEGIN_1 */
    result = Eigen::VectorXd::Zero(3);
    // midpoints of the edges of the reference triangle
    Eigen::MatrixXd midpoints_loc(2, 3);
    midpoints_loc << 0.5, 0.5, 0, 0, 0.5, 0.5;
    // obtain global midpoints
    auto midpoints_glob = geom->Global(midpoints_loc);
    // get determinants (scaling)
    const Eigen::VectorXd determinants(geom->IntegrationElement(midpoints_loc));
    // evaluate w(x) at the midpoints
    Eigen::VectorXd w(3);
    for (int i = 0; i < 3; i++) {
      w(i) = FESourceElemVecProvider::w_coeff_vec_(glob_indices[i]);
    }
    // add contributions to result
    for (int i = 0; i < 3; i++) {
      result(i) +=
          determinants(i) *
          (0.5 / (1. + std::pow(0.5 * w(i) + 0.5 * w((i + 1) % 3), 2)));
      result((i + 1) % 3) +=
          determinants(i) *
          (0.5 / (1. + std::pow(0.5 * w(i) + 0.5 * w((i + 1) % 3), 2)));
    }
    result *= 1. / 6.;
    /* SAM_LISTING_END_1 */
  } else if (ref_element.NumNodes() == 4) {
    /* SAM_LISTING_BEGIN_2 */
    result = Eigen::VectorXd::Zero(4);
    // midpoints of the edges of the unit square
    Eigen::MatrixXd midpoints_loc(2, 4);
    midpoints_loc << 0.5, 1, 0.5, 0, 0, 0.5, 1, 0.5;
    // obtain global midpoints
    auto midpoints_glob = geom->Global(midpoints_loc);
    // get determinants (scaling)
    const Eigen::VectorXd determinants(geom->IntegrationElement(midpoints_loc));
    // evaluate w(x) at the midpoints
    Eigen::VectorXd w(4);
    for (int i = 0; i < 4; i++) {
      w(i) = FESourceElemVecProvider::w_coeff_vec_(glob_indices[i]);
    }
    // add contributions to result
    for (int i = 0; i < 4; i++) {
      result(i) +=
          determinants(i) *
          (0.5 / (1. + std::pow(0.5 * w(i) + 0.5 * w((i + 1) % 4), 2)));
      result((i + 1) % 4) +=
          determinants(i) *
          (0.5 / (1. + std::pow(0.5 * w(i) + 0.5 * w((i + 1) % 4), 2)));
    }
    result *= 1. / 4.;
    /* SAM_LISTING_END_2 */
  } else {
    throw std::invalid_argument("received neither triangle nor quadrilateral");
  }
  /* SOLUTION_END */
  return result;
}
}  // namespace ParametricElementMatrices
