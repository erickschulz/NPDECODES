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
  /* SOLUTION_END */
  return result;
}
}  // namespace ParametricElementMatrices
