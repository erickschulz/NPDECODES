/**
 * @file
 * @brief NPDE homework SimpleLinearFiniteElements
 * @author Amélie Loher
 * @date 11/12/2019
 * @copyright Developed at ETH Zurich
 */

#include "simplelinearfiniteelements.h"

namespace SimpleLinearFiniteElements {

double getArea(const Eigen::Matrix<double, 2, 3>& triangle) {
  return 0.5 *
      ((triangle(0, 1) - triangle(0, 0)) * (triangle(1, 2) - triangle(1, 1)) -
       (triangle(0, 2) - triangle(0, 1)) * (triangle(1, 1) - triangle(1, 0)));
}

Eigen::Matrix<double, 2, 3> gradbarycoordinates(const Eigen::Matrix<double, 2, 3>& Vertices) {
  Eigen::Matrix3d X;

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
  double area = getArea(Vertices);
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
  double area = getArea(Vertices);
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
  double area = getArea(Vertices);
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
               const std::function<double(const Eigen::Vector2d&)> exact) {
  double l2error_squared = 0.0;
#if SOLUTION

  // loop over all triangles
  for (int i = 0; i < mesh.Elements.rows(); ++i) {
    Eigen::Matrix<double, 2, 3> triangle = mesh[i];

    // loop over all three vertices of the triangle
    Eigen::Vector3d error_at_vertices;
    for (int k = 0; k < 3; ++k) {
      error_at_vertices(k) = exact(triangle.col(k)) - uFEM(mesh.Elements(i, k));
    }

    // Add squared error per triangle
    l2error_squared += getArea(triangle) / 3.0 * error_at_vertices.squaredNorm();
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

  return std::sqrt(l2error_squared);
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
                const std::function<Eigen::Vector2d(const Eigen::Vector2d&)> exact) {
  double H1Serror_squared = 0.0;
#if SOLUTION

  // loop over all triangles
  for (int i = 0; i < mesh.Elements.rows(); ++i) {
    Eigen::Matrix<double, 2, 3> triangle = mesh[i];

    // loop over all three vertices of the triangle
    Eigen::Vector3d values_at_vertices;
    for (int k = 0; k < 3; ++k) {
      values_at_vertices(k) = uFEM[mesh.Elements(i, k)];
    }

    // gradient of FEM approximation (same for all 3 vertices!)
    Eigen::Vector2d gradientFEM = gradbarycoordinates(triangle) * values_at_vertices;

    // loop over all three vertices of the triangle
    Eigen::Vector3d error_at_vertices;
    for (int k = 0; k < 3; ++k) {
      Eigen::Vector2d gradient_exact = exact(triangle.col(k));
      error_at_vertices(k) = (gradientFEM - gradient_exact).squaredNorm();
    }

    // Add squared error per triangle
    H1Serror_squared += getArea(triangle) / 3.0 * error_at_vertices.sum();
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

  return std::sqrt(H1Serror_squared);
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
  auto uExact = [pi](const Eigen::Vector2d& x) {
    return std::cos(2 * pi * x(0)) * std::cos(2 * pi * x(1));
  };

  Eigen::VectorXd U;
  double l2error;
  double h1error;

#if SOLUTION
  // the gradient of uExact that can be easily analytically computed
  auto gradUExact = [pi](const Eigen::Vector2d& x) {
    Eigen::Vector2d gradient;
    gradient << -2 * pi * std::sin(2 * pi * x(0)) * std::cos(2 * pi * x(1)),
        -2 * pi * std::cos(2 * pi * x(0)) * std::sin(2 * pi * x(1));
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
