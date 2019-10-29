/* **********************************************************************
 Demo code for course "Numerical Methods for PDEs
 Section "Case Study: Triangular Linear FEM in Two Dimensions
 ********************************************************************** */

#include <mgl2/mgl.h>
#include "SimpleLinearFEM2D.h"

static mglGraph graph;

/**
 * @brief plots the mesh to an MathGL Plot
 * @param epsfile filepath, where plot should be stored as .eps file
 */
void plotMesh(const TriaMesh2D &mesh, const std::string &epsfile,
              int drawvertices = 0) {
  // get horizontal vertice data
  const Eigen::VectorXd &x = mesh.Coordinates.col(0);
  // get vertical  vertice data
  const Eigen::VectorXd &y = mesh.Coordinates.col(1);
  // get triangles indice data
  const Eigen::MatrixXi &T = mesh.Elements;
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
  // write .eps file to store graphics
  graph.WriteEPS(epsfile.c_str());
}

/**
 * @brief Creates a triangular surface plot based on the triangular mesh and the
 * provided values for the vertices
 * @param epsfile Filename where the plot should be stored as .eps file
 * @param values Vector containing a value for each vertex
 */
void plotSurf(const TriaMesh2D &mesh, const std::string &epsfile,
              const Eigen::VectorXd &values) {
  // get all x coordinates
  const Eigen::VectorXd &x = mesh.Coordinates.col(0);
  // get all y coordinates
  const Eigen::VectorXd &y = mesh.Coordinates.col(1);
  // get triangle indices
  const Eigen::MatrixXd &Tcopy =
      ((const Eigen::MatrixXd)mesh.Elements.cast<double>()).transpose();

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
