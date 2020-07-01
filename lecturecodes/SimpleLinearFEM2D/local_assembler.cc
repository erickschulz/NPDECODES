/* **********************************************************************
 Demo code for course "Numerical Methods for PDEs
 Section "Case Study: Triangular Linear FEM in Two Dimensions
 ********************************************************************** */

#include "local_assembler.h"

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 2, 3> gradbarycoordinates(const TriGeo_t& vertices) {
  Eigen::Matrix<double, 3, 3> X;
  // Argument \texttt{vertices} passes the vertex positions of the triangle
  // as the \textbf{rows} of a $3\times 2$-matrix, , see
  // \cref{cpp:getVtCoords}. The function returns the components of the
  // gradients as the \textbf{columns} of a $2\times 3$-matrix

  // Computation based on \eqref{eq:lambdalse}, solving for the
  // coefficients of the barycentric coordinate functions.
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  return X.inverse().block<2, 3>(1, 0);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::Matrix3d ElementMatrix_Lapl_LFE(const TriGeo_t& V) {
  // Argument \texttt{V} same as \texttt{vertices} in \cref{cpp:gradbarycoords}.
  // The function returns the $3\times 3$ element matrix as a fixed size
  // \eigen matrix.

  // Evaluate \eqref{eq:aK}, exploiting that the gradients are constant.
  // First compute the area of triangle by determinant formula
  double area = 0.5 * std::abs((V(0, 1) - V(0, 0)) * (V(1, 2) - V(1, 1)) -
                               (V(0, 2) - V(0, 1)) * (V(1, 1) - V(1, 0)));
  // Compute gradients of barycentric coordinate functions, see
  // \cref{cpp:gradbarycoords}
  Eigen::Matrix<double, 2, 3> X = gradbarycoordinates(V);
  // compute inner products of gradients through matrix multiplication
  return area * X.transpose() * X;
}
/* SAM_LISTING_END_2 */

Eigen::Matrix3d ElementMatrix_Mass_LFE(const TriGeo_t& vertices) {
  // compute area of triangle
  double area =
      0.5 *
      ((vertices(0, 1) - vertices(0, 0)) * (vertices(1, 2) - vertices(1, 1)) -
       (vertices(0, 2) - vertices(0, 1)) * (vertices(1, 1) - vertices(1, 0)));
  Eigen::Matrix3d X = Eigen::Matrix3d::Constant(area / 12.0);
  for (int i = 0; i < 3; i++) {
    X(i, i) *= 2.0;
  }
  return X;
}

Eigen::Vector3d localLoadLFE(const TriGeo_t& vertices,
                             const FHandle_t& FHandle) {
  // compute area of triangle, \emph{cf.} \cref{mc:ElementMatrix_LaplLFE}
  double area =
      0.5 *
      ((vertices(0, 1) - vertices(0, 0)) * (vertices(1, 2) - vertices(1, 1)) -
       (vertices(0, 2) - vertices(0, 1)) * (vertices(1, 1) - vertices(1, 0)));
  // evaluate source function for vertex location
  Eigen::Vector3d philoc = Eigen::Vector3d::Zero();
  for (int i = 0; i < 3; i++) {
    philoc(i) = FHandle(vertices.col(i));
  }
  // scale with $\frac{1}{3}\cdot$area of triangle
  philoc *= area / 3.0;
  return philoc;
}

Eigen::Matrix3d ElementMatrix_LaplMass_LFE(const TriGeo_t& Vertices) {
  return ElementMatrix_Lapl_LFE(Vertices) + ElementMatrix_Mass_LFE(Vertices);
}
