/* **********************************************************************
 Demo code for course "Numerical Methods for PDEs
 Section "Case Study: Triangular Linear FEM in Two Dimensions
 ********************************************************************** */

#include "SimpleLinearFEM2D.h"

Eigen::SparseMatrix<double> MatrixAssembler::Assemble(TriaMesh2D const &mesh) {
  // Get dimensions of the mesh
  int num_vertices = mesh.Coordinates.rows();
  int num_cells = mesh.Elements.rows();

  Triplet_t triplets;

  // Loop over elements and add local contributions.
  for (int i = 0; i < num_cells; i++) {
    // Get local$\to$global index mapping for current element
    Eigen::Vector3i element = mesh.Elements.row(i);
    TriGeo_t vertices;
    // Extract vertices of current element
    for (int j = 0; j < 3; j++) {
      vertices.col(j) = (mesh.Coordinates.row(element(j))).transpose();
    }
    // Compute element matrix
    Eigen::Matrix3d element_Matrix =
        MatrixAssembler::localMatrixHandle(vertices);
    // Add triplets from element matrix
    // Notice how we use the implicit local to global map given by the mesh data
    // structure to assign the parts of the element matrix to the correct
    // global indices of the galerkin matrix.
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        triplets.push_back({element(j), element(k), element_Matrix(j, k)});
      }
    }
  }
  // Build sparse matrix from triplets using \eigen's functions
  Eigen::SparseMatrix<double> A(num_vertices, num_vertices);
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();
  return A;
}

Eigen::SparseMatrix<double>
SlowMatrixAssembler::Assemble(TriaMesh2D const &mesh) {
  // Get dimensions of the mesh
  int num_vertices = mesh.Coordinates.rows();
  int num_cells = mesh.Elements.rows();

  Eigen::SparseMatrix<double> A(num_vertices, num_vertices);

  // Loop over elements and add local contributions.
  for (int i = 0; i < num_cells; i++) {
    // Get local$\to$global index mapping for current element
    Eigen::Vector3i element = mesh.Elements.row(i);
    TriGeo_t vertices;
    // Extract vertices of current element
    for (int j = 0; j < 3; j++) {
      vertices.col(j) = (mesh.Coordinates.row(element(j))).transpose();
    }
    // Compute element matrix
    Eigen::Matrix3d element_Matrix =
        SlowMatrixAssembler::localMatrixHandle(vertices);
    // Add triplets from element matrix
    // Notice how we use the implicit local to global map given by the mesh data
    // structure to assign the parts of the element matrix to the correct
    // global indices of the galerkin matrix.
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        A.coeffRef(element(j), element(k)) += element_Matrix(j, k);
      }
    }
  }
  // Build sparse matrix from triplets using \eigen's functions
  A.makeCompressed();
  return A;
}

/* Cell oriented assembly function, not using an auxiliary object */
Eigen::SparseMatrix<double>
assembleGalMatLFE(const TriaMesh2D &Mesh,
                  const LocalMatrixHandle_t getElementMatrix) {
  // Fetch the number of vertices
  int N = Mesh.Coordinates.rows();
  // Fetch the number of elements/cells, see \cref{par:trimesh2Ddata}
  int M = Mesh.Elements.rows();
  // Create empty sparse Galerkin matrix \cob{$\VA$}
  Eigen::SparseMatrix<double> A(N, N);
  // \com{Loop over elements} and ``distribute'' local contributions
  for (int i = 0; i < M; i++) {
    // Get local$\to$global index mapping for current element, \emph{cf.}
    // \eqref{eq:idxdef}
    Eigen::Vector3i dofhk = Mesh.Elements.row(i);
    TriGeo_t Vertices;
    // Extract vertices of current element, see \cref{par:trimesh2Ddata}
    for (int j = 0; j < 3; j++)
      Vertices.col(j) = Mesh.Coordinates.row(dofhk(j)).transpose();
    // Compute $3\times 3$ element matrix \cob{$\VA_{K}$}
    Eigen::Matrix3d Ak = getElementMatrix(Vertices);
    // Add local contribution to Galerkin matrix
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        A.coeffRef(dofhk(j), dofhk(k)) += Ak(j, k); // \Label[line]{asgl:1}
  }
  // Convert into CRS format, see \ncseref{sec:eigensparse}.
  A.makeCompressed();
  return A;
}
