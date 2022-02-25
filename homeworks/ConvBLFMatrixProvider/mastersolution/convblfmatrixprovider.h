#ifndef CDBLF_H_
#define CDBLF_H_
/**
 * @file
 * @brief demonstration of assembly of Galerkin linear system in LehrFEM++
 * assemble module; meant to provide sample codes for lecture document
 * @author Ralf Hiptmair
 * @date   April 2021
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

namespace cblfdemo {
// Auxiliary function computing the gradients of barycentric coordinate
// functions on a striangle and returning them in the columns of 2x3 matrix
std::pair<Eigen::Matrix<double, 2, 3>, double> getGradBaryCoords(
    const lf::mesh::Entity& tria);

/** @brief Compute element matrix for convection-diffusion bilinear form on
 * **triangles**
 */
template <typename MeshFunction>
class CDBLFElemMatProvider {
 public:
  // No default constructor
  CDBLFElemMatProvider() = delete;
  // Constructor takes a vectorfield argument
  explicit CDBLFElemMatProvider(MeshFunction av) : av_(std::move(av)) {}
  // Other constructors
  CDBLFElemMatProvider(const CDBLFElemMatProvider&) = default;
  CDBLFElemMatProvider(CDBLFElemMatProvider&&) noexcept = default;
  CDBLFElemMatProvider& operator=(const CDBLFElemMatProvider&) = delete;
  CDBLFElemMatProvider& operator=(CDBLFElemMatProvider&&) = delete;
  // The crucial interface methods for ENTITY_MATRIX_PROVIDER
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  Eigen::Matrix<double, 1, 3> Eval(const lf::mesh::Entity& tria);

 private:
  MeshFunction av_;
};

template <typename MeshFunction>
Eigen::Matrix<double, 1, 3> CDBLFElemMatProvider<MeshFunction>::Eval(
    const lf::mesh::Entity& tria) {
  // Throw error in case no triangular cell
  LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << tria.RefEl());
  // Fetch constant gradients of barycentric coordinate functions and area
  auto [grad_bary_coords, area] = getGradBaryCoords(tria);
  // Get values of vector field in midpoints of edges
  const Eigen::MatrixXd refqrnodes{
      (Eigen::MatrixXd(2, 3) << 0.5, 0.5, 0.0, 0.0, 0.5, 0.5).finished()};
  std::vector<Eigen::Vector2d> av_values{av_(tria, refqrnodes)};
  Eigen::Matrix<double, 1, 3> elmat{Eigen::Matrix<double, 1, 3>::Zero()};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      elmat(0, i) += (grad_bary_coords.col(i)).transpose() * av_values[j];
    }
  }
  return (area / 3) * elmat;
}

/** @brief Testing function for matrix provider for convection bilinear form
 *
 */
double testCDBLF(std::shared_ptr<lf::mesh::Mesh> mesh_p, Eigen::Vector2d a,
                 Eigen::Vector2d b);

}  // namespace cblfdemo

#endif
