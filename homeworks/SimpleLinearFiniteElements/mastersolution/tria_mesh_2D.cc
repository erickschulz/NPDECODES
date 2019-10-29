#include "simple_linear_finite_elements.h"

#include <iostream>

#include <boost/filesystem.hpp>

#include <mgl2/mgl.h>

namespace SimpleLinearFiniteElements {

static mglGraph graph;
/**
 * @brief read in mesh file into a costum data structure for this exercise
 * @param filename name of the mesh file
 */
TriaMesh2D::TriaMesh2D(const std::string &filename) {
  boost::filesystem::path here = __FILE__;
  auto filepath =
      (here.parent_path().parent_path() / ("meshes/" + filename)).string();

  std::ifstream mesh_file(filepath, std::ifstream::in);
  std::cout << "Load mesh " << filepath << std::endl;
  if (!mesh_file.good()) {
    throw std::runtime_error("Cannot open mesh file! File not found");
    return;
  }
  int nVertices;
  mesh_file >> nVertices;
  std::cout << nVertices << " Vertices" << std::endl;
  char keyword[1024];
  mesh_file.getline(keyword, 1024);
  if (!strcmp(keyword, "Vertices")) {
    throw std::runtime_error("Keyword 'Vertices' not found. Wrong file format");
    return;
  }
  Coordinates.resize(nVertices, 2);
  int nV = 0;
  while (nV < nVertices) {
    mesh_file >> Coordinates(nV, 0);
    mesh_file >> Coordinates(nV, 1);
    nV++;
  }
  int nElements;
  mesh_file >> nElements;
  mesh_file.getline(keyword, 1024);
  if (!strcmp(keyword, "Elements")) {
    throw std::runtime_error("Keyword 'Elements' not found. Wrong file format");
    return;
  }
  Elements.resize(nElements, 3);
  int nE = 0;
  while (nE < nElements) {
    mesh_file >> Elements(nE, 0);
    mesh_file >> Elements(nE, 1);
    mesh_file >> Elements(nE, 2);
    nE++;
  }
  mesh_file.close();
}

/**
 * @brief plots the mesh to an MathGL Plot
 * @param epsfile filepath, where plot should be stored as .eps file
 */
void TriaMesh2D::plotMesh(const std::string &epsfile,
                                int drawvertices) const {
  // get horizontal vertice data
  const Eigen::VectorXd &x = Coordinates.col(0);
  // get vertical  vertice data
  const Eigen::VectorXd &y = Coordinates.col(1);
  // get triangles indice data
  const Eigen::MatrixXi &T = Elements;
  const Eigen::MatrixXd &Tcopy =
      ((const Eigen::MatrixXd)T.cast<double>())
          .transpose();  // cast is needed, otherwise => errors
  // get data into mglData bins
  mglData xd(x.data(), x.size());
  mglData yd(y.data(), y.size());
  mglData Td(Tcopy.cols(), Tcopy.rows(), Tcopy.data());

  // give the plot an axis
  graph.SetRange('x', xd);
  graph.SetRange('y', yd);
  graph.Axis();
  // Plot each individual triangle
  graph.TriPlot(Td, xd, yd, "#b");
  if (drawvertices) {
    // Name the vertices from 1..n
    graph.Label(xd, yd, "%n");
    // Plot the vertices
    graph.Plot(xd, yd, " r*");
  }

  // write .eps file to storage
  graph.WriteEPS(epsfile.c_str());
}

/**
 * @brief Creates a triangular surface plot based on the triangular mesh and the
 * provided values for the vertices
 * @param epsfile Filename where the plot should be stored as .eps file
 * @param values Vector containing a value for each vertex
 */
void TriaMesh2D::plotSurf(const std::string &epsfile,
                                const Eigen::VectorXd &values) const {
  // get all x coordinates
  const Eigen::VectorXd &x = Coordinates.col(0);
  // get all y coordinates
  const Eigen::VectorXd &y = Coordinates.col(1);
  // get triangle indices
  const Eigen::MatrixXd &Tcopy =
      ((const Eigen::MatrixXd)Elements.cast<double>()).transpose();

  // prepare mglData bins
  mglData xd(x.data(), x.size());
  mglData yd(y.data(), y.size());
  mglData zd(values.size(), values.data());
  mglData Td(Tcopy.cols(), Tcopy.rows(), Tcopy.data());

  // plot
  mglGraph graph;
  // set appropriate axis ranges
  graph.SetRange('x', xd);
  graph.SetRange('y', yd);
  graph.SetRange('z', zd);
  // rotate for better visibility
  graph.Rotate(60, 45);
  graph.SetFontSizePT(8);
  graph.Grid("xyz", "h");
  graph.Axis();
  graph.Colorbar("");
  // plot and write to file
  graph.TriPlot(Td, xd, yd, zd);
  graph.WriteEPS(epsfile.c_str());
}

}  // namespace SimpleLinearFiniteElements
