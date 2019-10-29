/* **********************************************************************
 Demo code for course "Numerical Methods for PDEs
 Section "Case Study: Triangular Linear FEM in Two Dimensions
 ********************************************************************** */

#include "local_assembler.h"

Eigen::Matrix<double, 2, 3> gradbarycoordinates(const TriGeo_t& vertices) {
  Eigen::Matrix<double, 3, 3> X;

  // solve for the coefficients of the barycentric coordinate functions, see
  // \eqref{eq:lambdalse}
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  return X.inverse().block<2, 3>(1, 0);
}

Eigen::Matrix3d ElementMatrix_Lapl_LFE(const TriGeo_t& vertices) {
  // compute area of triangle
  double area =
      0.5 *
      ((vertices(0, 1) - vertices(0, 0)) * (vertices(1, 2) - vertices(1, 1)) -
       (vertices(0, 2) - vertices(0, 1)) * (vertices(1, 1) - vertices(1, 0)));
  // compute gradients of barycentric coordinate functions, see
  // \cref{mc:gradbarycoords}
  Eigen::Matrix<double, 2, 3> X = gradbarycoordinates(vertices);
  // compute inner products of gradients through matrix multiplication
  return area * X.transpose() * X;
}

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
