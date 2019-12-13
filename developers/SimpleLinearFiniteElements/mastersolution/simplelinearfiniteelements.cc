/**
 * @file
 * @brief NPDE homework SimpleLinearFiniteElements
 * @author Amélie Loher
 * @date 11/12/2019
 * @copyright Developed at ETH Zurich
 */

#include "simplelinearfiniteelements.h"

namespace SimpleLinearFiniteElements {

Eigen::Matrix<double, 2, 3> gradbarycoordinates(const Eigen::Matrix<double, 2, 3>& Vertices) {
  Eigen::Matrix<double, 3, 3> X;

  // solve for the coefficients of the barycentric coordinate functions, see
  // \eqref{eq:lambdalse}
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = Vertices.transpose();
  return X.inverse().block<2, 3>(1, 0);
}

/**
 *  @brief Computation of Element Matrix for the Laplacian
 */
Eigen::Matrix3d ElementMatrix_Lapl_LFE(const Eigen::Matrix<double, 2, 3>& Vertices) {
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
/**
 *  @brief Computation of full Galerkin Matrix
 */
Eigen::Matrix3d ElementMatrix_LaplMass_LFE(const Eigen::Matrix<double, 2, 3>& Vertices) {
  return ElementMatrix_Lapl_LFE(Vertices) + ElementMatrix_Mass_LFE(Vertices);
}

Eigen::Vector3d localLoadLFE(const Eigen::Matrix<double, 2, 3>& Vertices,
                             const std::function<double(const Eigen::Vector2d&)>& FHandle) {
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
 * 	@brief Computation of element mass matrix on planar triangle
 *  @param vertices 2x3 matrix of vertex coordinates
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix3d ElementMatrix_Mass_LFE(const Eigen::Matrix<double, 2, 3>& Vertices) {
  Eigen::Matrix3d result;
#if SOLUTION
  const Eigen::Vector2d a = (Vertices.col(1) - Vertices.col(0));
  const Eigen::Vector2d b = (Vertices.col(2) - Vertices.col(0));
  const double area = std::abs(a(0) * b(1) - a(1) * b(0)) / 2.;
  // initialize matrix
  result.setConstant(area / 12.);
  result += result.diagonal().asDiagonal();
#else
  //====================
  // Your code goes here
  //====================
#endif
  return result;
}
/* SAM_LISTING_END_1 */

/*double L2Norm_squared(const TriaMesh2D& mesh, const Eigen::VectorXd& u) {

  double l2norm_squared = 0.0;
  for (Eigen::Vector3d indices : mesh.Elements) {
    // vertices of triangle
    Eigen::Vector2d a = mesh.Vertices[indices(0)];
    Eigen::Vector2d b = mesh.Vertices[indices(1)];
    Eigen::Vector2d c = mesh.Vertices[indices(2)];

    // area of triangle
    double area = 0.5 * ((b(0) - a(0)) * (c(1) - a(1)) - (b(1) - a(1)) * (c(0) - a(0)));

    // values of u on vertices
    Eigen::Vector3d values_on_vertices(u[indices(0)], u[indices(1)], u[indices(2)]);

    // (approximate) L2-norm squared on single triangle
    l2norm_squared += area / 3.0 * values_on_vertices.squaredNorm();

    // TODO: Use e.g. Kahan summation algorithm to reduce cancellation
  }

  return l2norm_squared;
}*/

/**
 * @brief L2Error Computes the L2 error between the approximate solution and
 *                the exact solution
 * @param mesh the mesh to use
 * @param uFEM the solution approximated through FEM
 * @param exact the exact solution
 * @return the L2 difference
 */
/* SAM_LISTING_BEGIN_2 */
double L2Error(const TriaMesh2D& mesh, const Eigen::VectorXd& uFEM,
               const std::function<double(double, double)> exact) {
  double error;
  const auto& triangles = mesh.Elements;
  const auto& vertices = mesh.Coordinates;

#if SOLUTION
  for (size_t triangle = 0; triangle < triangles.rows(); ++triangle) {
    // global vertex numbers
    const auto& indexSet = triangles.row(triangle);

    // get the 3 vertices of the triangle
    const auto& a = vertices.row(indexSet(0));
    const auto& b = vertices.row(indexSet(1));
    const auto& c = vertices.row(indexSet(2));

    // compute area of triangle by determinant formula
    const double area =
        0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));

    // compute integral over triangle by quadrature
    double localError = 0;
    for (size_t vertex = 0; vertex < 3; ++vertex) {
      const auto& v = vertices.row(indexSet(vertex));
      const double fe_val = uFEM[indexSet(vertex)];
      const double exact_val = exact(v[0], v[1]);
      localError += std::pow(fe_val - exact_val, 2);
   } 
    localError *= (area / 3.0);
    // Sum up local error contributions
    error += localError;
  }  // end of main loop over triangles
#else
  //====================
  // Your code goes here
  //====================
#endif
  return std::sqrt(error);
} 
/* SAM_LISTING_END_2 */

/**
 * @brief H1Serror Computes the H^1 error between the approximate solution and
 *                the exact solution
 * @param mesh the mesh to use
 * @param uFEM the solution approximated through FEM
 * @param exact the exact gradient of the solution
 * @return the H^1 difference
 *
 * @note This implementation seems to be flawed!
 */
/* SAM_LISTING_BEGIN_3 */
double H1Serror(const TriaMesh2D& mesh, const Eigen::VectorXd& uFEM,
                const std::function<Eigen::Vector2d(double, double)> exact) {
  double error = 0;
#if SOLUTION
  const auto& triangles = mesh.Elements;
  const auto& vertices = mesh.Coordinates;

  for (size_t triangle = 0; triangle < triangles.rows(); ++triangle) {
    // global vertex numbers
    const auto& indexSet = triangles.row(triangle);

    // get the 3 vertices of the triangle
    const auto& a = vertices.row(indexSet(0));
    const auto& b = vertices.row(indexSet(1));
    const auto& c = vertices.row(indexSet(2));

    // compute area of triangle
    double area =
        0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));

    double localError = 0;
    Eigen::Matrix<double, 2, 3> element;
    // Extract vertices of current element, see \lref{par:trimesh2Ddata}
    for (int j = 0; j < 3; j++)
      element.col(j) = mesh.Coordinates.row(indexSet(j)).transpose();
    
	auto gradBarycentric = gradbarycoordinates(element);

    // integrate the exact solution over the element by quadrature
    // the finite element solution need not be aproximated as it is constant
    Eigen::Vector2d gradientFEM;
    gradientFEM.setZero();
    // compute constant gradient over the triangle
    for (size_t vertex = 0; vertex < 3; ++vertex) {
      gradientFEM += uFEM[indexSet(vertex)] * gradBarycentric.col(vertex);
    }
    
	for (size_t vertex = 0; vertex < 3; ++vertex) {
      const auto& v = vertices.row(indexSet(vertex));
      Eigen::VectorXd gradientExact = exact(v[0], v[1]);
      for (int j = 0; j < 2; ++j)
        localError += std::pow(gradientFEM(j) - gradientExact(j), 2);
    }
    localError *= (area / 3.);
    // Sum up local error contributions
    error += localError;
  }  // end of main loop over triangles
#else
  //====================
  // Your code goes here
  //====================
#endif
  return std::sqrt(error);
}
/* SAM_LISTING_END_3 */

/**
 * @brief assemLoad_LFE Assembles the Load Vector
 * @param mesh the mesh to use
 * @param getElementVector 
 * @param FHandle function handle for f
 * @return assembled load vector
 */
Eigen::VectorXd assemLoad_LFE(const TriaMesh2D& Mesh,
                              const std::function<Eigen::Vector3d(const Eigen::Matrix<double, 2, 3>&,
							  std::function<double(const Eigen::Vector2d&)>)>& getElementVector,
                              const std::function<double(const Eigen::Vector2d&)>& FHandle) {
  // obtain the number of vertices
  int N = Mesh.Coordinates.rows();
  // obtain the number of elements/cells
  int M = Mesh.Elements.rows();
  Eigen::VectorXd phi = Eigen::VectorXd::Zero(N);
  // loop over elements and add local contributions
  for (int i = 0; i < M; i++) {
    // get local$\to$global index mapping for current element
    Eigen::Vector3i element = Mesh.Elements.row(i);
    Eigen::Matrix<double, 2, 3> Vertices;
    // extract vertices of current element
    for (int j = 0; j < 3; j++) {
      Vertices.col(j) = Mesh.Coordinates.row(element(j)).transpose();
    }
    // compute element right hand side vector
    Eigen::Vector3d philoc = getElementVector(Vertices, FHandle);
    // add contributions to global load vector
    for (int j = 0; j < 3; j++) {
      phi(element(j)) += philoc(j);
    }
  }
  return phi;
}

/**
 * @brief assemLoad_LFE Assembles the Load Vector
 * @param mesh the mesh to use
 * @param getElementMatrix Element Matrix  
 * @return Galerkin Matrix
 */
Eigen::SparseMatrix<double> GalerkinAssembly(
    const TriaMesh2D& Mesh, 
	const std::function<Eigen::Matrix3d(const Eigen::Matrix<double, 2, 3>&)>& getElementMatrix) {
  
  // obtain the number of vertices
  int N = Mesh.Coordinates.rows();
  // obtain the number of elements/cells
  int M = Mesh.Elements.rows();
  std::vector<Eigen::Triplet<double> > triplets;
  // loop over elements and add local contributions
  for (int i = 0; i < M; i++) {
    // get local$\to$global index mapping for current element, \emph{cf.}
    // \lref{eq:idxdef}
    Eigen::Vector3i element = Mesh.Elements.row(i);
    Eigen::Matrix<double, 2, 3> Vertices;
    // extract vertices of current element
    for (int j = 0; j < 3; j++) {
      Vertices.col(j) = Mesh.Coordinates.row(element(j)).transpose();
    }
    // compute element contributions
    Eigen::Matrix3d Ak = getElementMatrix(Vertices);
    // build triplets from contributions
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        triplets.push_back({element(j), element(k), Ak(j, k)});
      }
    }
  }
  // build sparse matrix from triplets
  Eigen::SparseMatrix<double> A(N, N);
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();
  return A;
}


/**
 * @brief solves system and prints H1-semierror, L2 error, the mesh and a
 * surface plot
 * @param mesh: discretisation of the computational domain
 */
/* SAM_LISTING_BEGIN_4 */
std::tuple<Eigen::VectorXd, double, double> solve(const SimpleLinearFiniteElements::TriaMesh2D &mesh) {
  
  const double pi = 3.1415926535897;

  // define the source function f
  std::function<double(const Eigen::Vector2d&)> f = [pi](const Eigen::Vector2d &x) {
	return (8.0 * pi * pi + 1) * std::cos(2 * pi * x(0)) *
           std::cos(2 * pi * x(1));
  };

  // the exact solution of the linear variational problem
  auto uExact = [pi](double x, double y) {
    return std::cos(2 * pi * x) * std::cos(2 * pi * y);
  };

  Eigen::VectorXd U;
  double l2error;
  double h1error;

#if SOLUTION
  // the gradient of uExact that can be easily analytically computed
  auto gradUExact = [pi](double x, double y) {
    Eigen::Vector2d gradient;
    gradient << -2 * pi * std::sin(2 * pi * x) * std::cos(2 * pi * y),
        -2 * pi * std::cos(2 * pi * x) * std::sin(2 * pi * y);
    return gradient;
  };
  // Compute the galerkin matrix A and load vector L
  Eigen::SparseMatrix<double> A = SimpleLinearFiniteElements::GalerkinAssembly(
      mesh, SimpleLinearFiniteElements::ElementMatrix_LaplMass_LFE);

  Eigen::VectorXd L = SimpleLinearFiniteElements::assemLoad_LFE(
      mesh, SimpleLinearFiniteElements::localLoadLFE, f);

  // solve the LSE using the sparse LU solver of Eigen
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  U = solver.solve(L);

  // compute both the L2 and the H1s error
  l2error = SimpleLinearFiniteElements::L2Error(mesh, U, uExact);
  h1error = SimpleLinearFiniteElements::H1Serror(mesh, U, gradUExact);
#else
  //====================
  // Your code goes here
  // Assigning some dummy values
  U = Eigen::VectorXd::Zero(mesh.Coordinates.rows());
  l2error = 1.0;
  h1error = 1.0;
  //====================
#endif

  return std::make_tuple(U, l2error, h1error);
}
/* SAM_LISTING_END_4 */

} // namespace SimpleLinearFiniteElements 
