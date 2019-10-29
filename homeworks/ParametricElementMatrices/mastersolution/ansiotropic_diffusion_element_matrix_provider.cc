/**
 * @file ansiotropic_diffusion_element_matrix_provider.cc
 * @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans
 * @date 13/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "ansiotropic_diffusion_element_matrix_provider.h"

namespace ParametricElementMatrices {

/**
 * @brief compute element matrix for \int_{\Omega} grad(u(x)) (1 + d(x)d(x)')
 * grad(v(x)) dx
 * @param cell edge on the boundary
 */
AnisotropicDiffusionElementMatrixProvider::ElemMat
AnisotropicDiffusionElementMatrixProvider::Eval(const lf::mesh::Entity &cell) {
  auto geom = cell.Geometry();
  auto ref_element = geom->RefEl();
  Eigen::MatrixXd result;
  /* SOLUTION_BEGIN */

  // distinguish between triangles and quadrilaterals
  if (ref_element.NumNodes() == 3) {
    /* SAM_LISTING_BEGIN_1 */
    result = Eigen::MatrixXd::Zero(3, 3);
    // midpoints of the edges of the reference triangle
    Eigen::Matrix<double, 2, 3> midpoints_loc(2, 3);
    midpoints_loc << 0.5, 0.5, 0, 0, 0.5, 0.5;
    // obtain global midpoints
    auto midpoints_glob = geom->Global(midpoints_loc);

    // get jacobian and determinants (scaling)
    const Eigen::MatrixXd JinvT(geom->JacobianInverseGramian(midpoints_loc));
    const Eigen::VectorXd determinants(geom->IntegrationElement(midpoints_loc));
    // local gradients
    Eigen::Matrix<double, 2, 3> grads_loc(2, 3);
    grads_loc << -1, 1, 0, -1, 0, 1;
    // loop over quadrature points
    for (int i = 0; i < 3; i++) {
      // prepare extra factor matrix
      Eigen::Vector2d d_x = AnisotropicDiffusionElementMatrixProvider::Vf_d_(
          midpoints_glob.col(i));
      Eigen::Matrix2d factor_matrix =
          Eigen::Matrix2d::Identity(2, 2) + d_x * d_x.transpose();
      // compute global gradients
      Eigen::MatrixXd grad_glob(2, 3);
      grad_glob = JinvT.block(0, 2 * i, 2, 2) * grads_loc;
      // add contributions to element matrix
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          double tmp =
              grad_glob.col(j).transpose() * factor_matrix * grad_glob.col(k);
          result(j, k) += determinants(i) * tmp;
        }
      }
    }
    // multiply with area divided by 6 (3 quad points and RefEl has area 0.5)
    result *= 1. / 6.;
    /* SAM_LISTING_END_1 */
  } else if (ref_element.NumNodes() == 4) {
    /* SAM_LISTING_BEGIN_2 */
    result = Eigen::MatrixXd::Zero(4, 4);
    // midpoints of the edges of the reference triangle
    Eigen::Matrix<double, 2, 4> midpoints_loc(2, 4);
    midpoints_loc << 0.5, 1, 0.5, 0, 0, 0.5, 1, 0.5;
    // obtain global midpoints
    auto midpoints_glob = geom->Global(midpoints_loc);
    // get jacobian and determinants (scaling)
    const Eigen::MatrixXd JinvT(geom->JacobianInverseGramian(midpoints_loc));
    const Eigen::VectorXd determinants(geom->IntegrationElement(midpoints_loc));
    // analytic gradients for all for basis functions
    auto grad_f = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 4> {
      Eigen::Matrix<double, 2, 4> result;
      result(0, 0) = x(1) - 1;
      result(1, 0) = x(0) - 1;
      result(0, 1) = 1 - x(1);
      result(1, 1) = -x(0);
      result(0, 2) = x(1);
      result(1, 2) = x(0);
      result(0, 3) = -x(1);
      result(1, 3) = 1 - x(0);
      return result;
    };
    // Loop over quadrature points
    for (int i = 0; i < 4; i++) {
      // prepare extra factor matrix
      Eigen::Vector2d d_x = AnisotropicDiffusionElementMatrixProvider::Vf_d_(
          midpoints_glob.col(i));
      Eigen::Matrix2d factor_matrix =
          Eigen::Matrix2d::Identity(2, 2) + d_x * d_x.transpose();
      // compute global gradients
      Eigen::Matrix<double, 2, 4> grad_glob{JinvT.block(0, 2 * i, 2, 2) *
                                            grad_f(midpoints_loc.col(i))};
      // add contributions to element matrix
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          double tmp =
              (grad_glob.col(j).transpose() * factor_matrix * grad_glob.col(k));
          result(j, k) += determinants(i) * tmp;
        }
      }
    }
    // multiply with area divided by 4 (number of quad points)
    result *= 1. / 4.;
    /* SAM_LISTING_END_2 */
  } else {
    throw std::invalid_argument("received neither triangle nor quadrilateral");
  }
  /* SOLUTION_END */
  return result;
}

/**
 * @brief constructor of AnisotropicDiffusionElementMatrixProvider for
 * \int_{\Omega} grad(u(x)) (1 + d(x)d(x)') grad(v(x)) dx
 * @param Vf_d vectorfield d(x)
 */
AnisotropicDiffusionElementMatrixProvider::
    AnisotropicDiffusionElementMatrixProvider(
        AnisotropicDiffusionElementMatrixProvider::vectorfield_t Vf_d) {
  AnisotropicDiffusionElementMatrixProvider::Vf_d_ = Vf_d;
}
}  // namespace ParametricElementMatrices
