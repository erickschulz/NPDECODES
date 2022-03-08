/* **********************************************************************
   Demo code for course "Numerical Methods for PDEs
   Section "Case Study: Triangular Linear FEM in Two Dimensions
   ********************************************************************** */

#include "SimpleLinearFEM2D.h"

Eigen::VectorXd VectorAssembler::Assemble(TriaMesh2D const& mesh) {
  // obtain the number of vertices
  int num_vertices = mesh._nodecoords.rows();
  // obtain the number of cells
  int num_cells = mesh._elements.rows();

  Eigen::VectorXd phi = Eigen::VectorXd::Zero(num_vertices);

  for (int i = 0; i < num_cells; i++) {
    // get local to global map
    Eigen::Vector3i element = mesh._elements.row(i);
    // get vertices of current triangle
    TriGeo_t vertices;
    for (int j = 0; j < 3; j++) {
      vertices.col(j) = (mesh._nodecoords.row(element(j))).transpose();
    }
    // compute element right hand side vector
    Eigen::Vector3d phi_loc = VectorAssembler::localVectorHandle(
        vertices, VectorAssembler::sourceFunction);
    for (int j = 0; j < 3; j++) phi(element(j)) += phi_loc(j);
  }
  return phi;
}
