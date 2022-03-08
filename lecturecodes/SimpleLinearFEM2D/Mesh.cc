/* **********************************************************************
   Demo code for course "Numerical Methods for PDEs
   Section "Case Study: Triangular Linear FEM in Two Dimensions
   ********************************************************************** */
#ifndef LINEMAX
#define LINEMAX 1024
#endif

#include "SimpleLinearFEM2D.h"

/**
 * Constructs TriaMesh2D.
 *
 * @param filename  location of the desired mesh file on disk
 */
TriaMesh2D::TriaMesh2D(std::string filename) {
  // Open mesh file
  ifstream mesh_file(filename, ifstream::in);
  cout << "Loading mesh from file " << filename << endl;
  if (!mesh_file.good()) {
    throw runtime_error("Cannot open mesh file!");
    return;
  }
  int num_vertices;
  mesh_file >> num_vertices;  // Read no. of nodes
  cout << num_vertices << " Vertices" << endl;
  char keyword[LINEMAX];
  mesh_file.getline(keyword, LINEMAX);
  if (!strcmp(keyword, "Vertices")) {
    throw runtime_error("Keyword 'Vertices' not found.");
    return;
  }
  // Read coordinates of nodes
  _nodecoords.resize(num_vertices, 2);
  for (int i = 0; i < num_vertices; i++) {
    mesh_file >> _nodecoords(i, 0);
    mesh_file >> _nodecoords(i, 1);
  }
  int num_elements;
  mesh_file >> num_elements;  // Read no. of cells
  mesh_file.getline(keyword, LINEMAX);
  if (!strcmp(keyword, "Elements")) {
    throw runtime_error("Keyword 'Elements' not found.");
    return;
  }
  // Read node indices of triangles
  _elements.resize(num_elements, 3);
  for (int i = 0; i < num_elements; i++) {
    mesh_file >> _elements(i, 0);
    mesh_file >> _elements(i, 1);
    mesh_file >> _elements(i, 2);
  }
  mesh_file.close();
}

/**
 *Retrieve coordinates of vertices of a triangles as rows of a 3x2
 *matrix
 *
 *@param cell_index    Index of a cell C
 *@return              Coordinates of vertices of C
 */
/* SAM_LISTING_BEGIN_7 */
TriGeo_t TriaMesh2D::getVtCoords(size_t cell_index) const {
  // Check whether valid cell index (starting from zero!)
  assert(cell_index < _elements.rows());
  // Obtain numbers of vertices of triangle i
  Eigen::RowVector3i idx = _elements.row(cell_index);
  // Bild matrix of vertex coordinates
  Eigen::Matrix<double, 3, 2> vtc;
  vtc << _nodecoords.row(idx[0]), _nodecoords.row(idx[1]),
      _nodecoords.row(idx[2]);
  return vtc.transpose();
}
/* SAM_LISTING_END_7 */
