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
