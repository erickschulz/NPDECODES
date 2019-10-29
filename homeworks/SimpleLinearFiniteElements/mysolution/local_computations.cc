#include "simple_linear_finite_elements.h"

#include <Eigen/Dense>

namespace SimpleLinearFiniteElements {

/**
 * @brief Computation of local load vector for linear Lagrangian finite elements
 * and linear form $\int fv dx$. Local approximation of
 * integral by 2D trapezoidal rule
 * @param Vertices Vertex positions of the triangle as the rows of a $3\times
 * 2$-matrix
 * @param FHandle A handle to a function that provides f(x) for the
 * source function f and for any point x in Omega
 * @return local load vector
 */
Eigen::Vector3d localLoadLFE(const TriGeo_t& Vertices,
                             const FHandle_t& FHandle) {
  // compute area of triangle, \emph{cf.} \cref{mc:ElementMatrix_LaplLFE}
  double area =
      0.5 *
      ((Vertices(0, 1) - Vertices(0, 0)) * (Vertices(1, 2) - Vertices(1, 1)) -
       (Vertices(0, 2) - Vertices(0, 1)) * (Vertices(1, 1) - Vertices(1, 0)));
  // evaluate source function for vertex location
  Eigen::Vector3d philoc = Eigen::Vector3d::Zero();
  for (int i = 0; i < 3; i++) {
    philoc(i) = FHandle(Vertices.col(i));
  }
  // scale with $\frac{1}{3}\cdot$area of triangle
  philoc *= area / 3.0;
  return philoc;
}

/**
 * @brief Function computing the gradients of barycentric coordinate functions
 * on a triangle.
 * @param Vertices Vertex positions of the triangle as the rows of a $3\times
 * 2$-matrix.
 * @return The components of the gradients as the columns of a $2\times
 * 3$-matrix.
 */
Eigen::Matrix<double, 2, 3> gradbarycoordinates(const TriGeo_t& Vertices) {
  Eigen::Matrix<double, 3, 3> X;

  // solve for the coefficients of the barycentric coordinate functions, see
  // \eqref{eq:lambdalse}
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = Vertices.transpose();
  return X.inverse().block<2, 3>(1, 0);
}

/**
 * @brief Function computing the stiffness element matrix for piecewise linear
 * Lagrangian finite elements on a triangle
 * @param Vertices Vertex positions of the triangle as the rows of a $3\times
 * 2$-matrix
 * @return The stiffness element matrix
 */
Eigen::Matrix3d ElementMatrix_Lapl_LFE(const TriGeo_t& Vertices) {
  // compute area of triangle
  double area =
      0.5 *
      ((Vertices(0, 1) - Vertices(0, 0)) * (Vertices(1, 2) - Vertices(1, 1)) -
       (Vertices(0, 2) - Vertices(0, 1)) * (Vertices(1, 1) - Vertices(1, 0)));
  // compute gradients of barycentric coordinate functions
  Eigen::Matrix<double, 2, 3> X = gradbarycoordinates(Vertices);
  // compute inner products of gradients through matrix multiplication
  return area * X.transpose() * X;
}

Eigen::Matrix3d ElementMatrix_LaplMass_LFE(const TriGeo_t& Vertices) {
  return ElementMatrix_Lapl_LFE(Vertices) + ElementMatrix_Mass_LFE(Vertices);
}

}  // namespace SimpleLinearFiniteElements
