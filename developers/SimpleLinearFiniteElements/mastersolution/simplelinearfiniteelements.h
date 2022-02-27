/**
 * @file
 * @brief NPDE homework SimpleLinearFiniteElements
 * @author Amélie Loher
 * @date 11/12/2019
 * @copyright Developed at ETH Zurich
 */

#ifndef SIMPLELINEARFINITEELEMENTS_H_
#define SIMPLELINEARFINITEELEMENTS_H_

#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <cmath>
#include <functional>
#include <string>
#include <tuple>

#include "tria_mesh_2D.h"

namespace SimpleLinearFiniteElements {

double getArea(const Eigen::Matrix<double, 2, 3> &triangle);

Eigen::Matrix<double, 2, 3> gradbarycoordinates(
    const Eigen::Matrix<double, 2, 3> &triangle);

Eigen::Matrix3d ElementMatrix_Mass_LFE(
    const Eigen::Matrix<double, 2, 3> &triangle);

Eigen::Matrix3d ElementMatrix_Lapl_LFE(
    const Eigen::Matrix<double, 2, 3> &triangle);

Eigen::Matrix3d ElementMatrix_LaplMass_LFE(
    const Eigen::Matrix<double, 2, 3> &triangle);

/**
 * @brief L2Error Computes the L2 error between the approximate solution and
 *                the exact solution
 * @tparam FUNCTOR to be compatible with
 *                 std::function<double(const Eigen::Vector2d &)>
 * @param mesh the mesh to use
 * @param uFEM the solution approximated through FEM
 * @param exact functor object providing the exact solution
 * @return the L2 difference
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
double L2Error(const TriaMesh2D &mesh, const Eigen::VectorXd &uFEM,
               FUNCTOR &&exact) {
  double l2error_squared = 0.0;
#if SOLUTION
  // loop over all triangles
  for (int i = 0; i < mesh.elements.rows(); ++i) {
    Eigen::Matrix<double, 2, 3> triangle = mesh[i];

    // loop over all three vertices of the triangle
    Eigen::Vector3d error_at_vertices;
    for (int k = 0; k < 3; ++k) {
      error_at_vertices(k) = exact(triangle.col(k)) - uFEM(mesh.elements(i, k));
    }
    // Add squared error per triangle
    l2error_squared +=
        getArea(triangle) / 3.0 * error_at_vertices.squaredNorm();
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return std::sqrt(l2error_squared);
}
/* SAM_LISTING_END_2 */

double L2Error_old(const SimpleLinearFiniteElements::TriaMesh2D &mesh,
                   const Eigen::VectorXd &uFEM,
                   const std::function<double(const Eigen::Vector2d &)> exact);

double H1Serror(
    const SimpleLinearFiniteElements::TriaMesh2D &mesh,
    const Eigen::VectorXd &uFEM,
    const std::function<Eigen::Vector2d(const Eigen::Vector2d &)> exact);

Eigen::SparseMatrix<double> GalerkinAssembly(
    const SimpleLinearFiniteElements::TriaMesh2D &mesh,
    const std::function<Eigen::Matrix3d(const Eigen::Matrix<double, 2, 3> &)>
        &getElementMatrix);

Eigen::VectorXd assemLoad_LFE(
    const SimpleLinearFiniteElements::TriaMesh2D &mesh,
    const std::function<double(const Eigen::Vector2d &)> &f);

std::tuple<Eigen::VectorXd, double, double> Solve(
    const SimpleLinearFiniteElements::TriaMesh2D &mesh);

}  // namespace SimpleLinearFiniteElements

#endif  // SIMPLELINEARFINITEELEMENTS_H_
