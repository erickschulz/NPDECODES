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
  // In order to avoid cancellation, use Taylor polynomial approximation of
  // the denominator for small arguments
  if (std::abs(tau) < 1e-10) {
    res = 1.0;
  } else if (std::abs(tau) < 1e-3) {
    res = 1.0 / (1.0 + (0.5 + 1.0 / 6.0 * tau) * tau);
  } else {
    // No cancellation for large arguments: use formula defining the Bernoulli
    // function
    res = tau / (std::exp(tau) - 1.0);
  }
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
  // Get the DOF handler for the FE space
  const lf::assemble::DofHandler& dofh = fe_space.LocGlobMap();
  // compute beta(e) for all edges of the mesh
  for (const lf::mesh::Entity* edge : fe_space.Mesh()->Entities(1)) {
    // Compute the DOF indices of the vertices of the edge
    auto idxs = dofh.GlobalDofIndices(*edge);
    (*beta_p)(*edge) =
        std::exp(mu(idxs[1])) * Bernoulli(mu(idxs[1]) - mu(idxs[0]));
  }

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
  // Get the values of beta on the edges of the triangle.
  // by the Lehrfem++ numbering convention
  // b = [beta(e_0), beta(e_1), beta(e_2)]' = [\beta_{1,2}, \beta_{2,3},
  // \beta_{1,3}]'
  Eigen::Vector3d b = beta_loc(cell);

  // evaluate the element matrix using the formula in subproblem h)
  result << AK(0, 1) * b(0) + AK(0, 2) * b(2), -AK(0, 1) * b(0),
      -AK(0, 2) * b(2), -AK(0, 1) * b(0), AK(0, 1) * b(0) + AK(1, 2) * b(1),
      -AK(1, 2) * b(1), -AK(0, 2) * b(2), -AK(1, 2) * b(1),
      AK(0, 2) * b(2) + AK(1, 2) * b(1);

  Eigen::Vector3d mu_exp = (-mu_loc(cell)).array().exp();
  result *= mu_exp.asDiagonal();
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
