/**
 * @file expfittedupwind.cc
 * @brief NPDE homework ExpFittedUpwind
 * @author Amélie Loher, Philippe Peter
 * @date 07.01.2021
 * @copyright Developed at ETH Zurich
 */

#include "expfittedupwind.h"

#include <lf/base/base.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

#include <Eigen/Core>
#include <cmath>
#include <memory>
#include <vector>

namespace ExpFittedUpwind {

/**
 * @brief Computes the Bernoulli function B(tau)
 **/
/* SAM_LISTING_BEGIN_1 */
double Bernoulli(double tau) {
  double res = 0;
  //====================
  // Your code goes here
  //====================
  return res;
}
/* SAM_LISTING_END_1 */

/**
 * @brief computes the quantities \beta(e) for all the edges e of a mesh
 * @param mesh_p underlying mesh
 * @param mu vector of nodal values of a potential Psi
 * @return  Mesh Data set containing the quantities \beta(e)
 */
/* SAM_LISTING_BEGIN_2 */
// REVISE: Does not match specification
std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<double>> CompBeta(
    const lf::uscalfe::FeSpaceLagrangeO1<double>& fe_space,
    const Eigen::VectorXd& mu) {
  // data set over all edges of the mesh.
  auto beta_p = lf::mesh::utils::make_CodimMeshDataSet(fe_space.Mesh(), 1, 1.0);
  //====================
  // Your code goes here
  //====================

  return beta_p;
}
/* SAM_LISTING_END_2 */

/**
 * @brief actual computation of the element matrix
 * @param cell reference to the triangle for which the matrix is evaluated
 * @return 3x3 dense matrix containg the element matrix
 */
/* SAM_LISTING_BEGIN_3 */
Eigen::Matrix3d ExpFittedEMP::Eval(const lf::mesh::Entity& cell) {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Only 2D triangles are supported.");

  // Evaluate the element matrix A_K
  Eigen::Matrix3d AK = laplace_provider_.Eval(cell).block<3, 3>(0, 0);
  Eigen::Matrix3d result;
  //====================
  // Your code goes here
  //====================
  return std::move(result);
}
/* SAM_LISTING_END_3 */

/**
 * @brief returns the quanties beta(e) for  the
 * three edges e_0, e_1 and e_2 of a triangle.
 * @param cell reference to the triangle for which the quantities are needed
 * @return vector  [beta(e_0),beta(e_1),beta(e_2)]'
 **/
Eigen::Vector3d ExpFittedEMP::beta_loc(const lf::mesh::Entity& cell) {
  Eigen::Vector3d b;
  auto edges = cell.SubEntities(1);
  for (int i = 0; i < 3; ++i) {
    b(i) = (*beta_)(*(edges[i]));
  }
  return b;
}

/** @brief returns the nodal values of the potential Psi for the
 * three vertices a_1, a_2 and a_3 of a triangle
 * @param cell reference to the triangle for which the quantities are needed
 * @return vector [Psi(a_1), Psi(a_2), Psi(a_3)]'
 **/
Eigen::Vector3d ExpFittedEMP::mu_loc(const lf::mesh::Entity& cell) {
  Eigen::Vector3d m;
  auto mesh_p = fe_space_->Mesh();
  auto vertices = cell.SubEntities(2);
  for (int i = 0; i < 3; ++i) {
    int index = mesh_p->Index(*(vertices[i]));
    m(i) = mu_(index);
  }
  return m;
}

} /* namespace ExpFittedUpwind */
