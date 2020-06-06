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
TriaMesh2D::TriaMesh2D(const string &filename) {
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
  Coordinates.resize(num_vertices, 2);
  for (int i = 0; i < num_vertices; i++) {
    mesh_file >> Coordinates(i, 0);
    mesh_file >> Coordinates(i, 1);
  }
  int num_elements;
  mesh_file >> num_elements;  // Read no. of cells
  mesh_file.getline(keyword, LINEMAX);
  if (!strcmp(keyword, "Elements")) {
    throw runtime_error("Keyword 'Elements' not found.");
    return;
  }
  // Read node indices of triangles
  Elements.resize(num_elements, 3);
  for (int i = 0; i < num_elements; i++) {
    mesh_file >> Elements(i, 0);
    mesh_file >> Elements(i, 1);
    mesh_file >> Elements(i, 2);
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
  assert(cell_index < Elements.rows());
  // Obtain numbers of vertices of triangle i
  Eigen::RowVector3i idx = Elements.row(cell_index);
  // Bild matrix of vertex coordinates
  Eigen::Matrix<double, 3, 2> vtc;
  vtc << Coordinates.row(idx[0]), Coordinates.row(idx[1]),
    Coordinates.row(idx[2]);
  return vtc.transpose();
}
/* SAM_LISTING_END_7 */
