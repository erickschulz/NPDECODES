
/**
 * @file ansiotropic_diffusion_element_matrix_provider.h
 * @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans
 * @date 13/03/2019
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

class AnisotropicDiffusionElementMatrixProvider {
 public:
  using elem_mat_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  using ElemMat = const elem_mat_t;
  using vectorfield_t = std::function<Eigen::Vector2d(Eigen::Vector2d)>;
  AnisotropicDiffusionElementMatrixProvider(vectorfield_t Vf_d);
  bool isActive(const lf::mesh::Entity &) { return true; }
  ElemMat Eval(const lf::mesh::Entity &cell);

 private:
  vectorfield_t Vf_d_;
};

}  // namespace ParametricElementMatrices
