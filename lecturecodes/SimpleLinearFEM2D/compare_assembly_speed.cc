/* **********************************************************************
 Demo code for course "Numerical Methods for PDEs
 Section "Case Study: Triangular Linear FEM in Two Dimensions
 ********************************************************************** */

#include <mgl2/mgl.h>
#include "SimpleLinearFEM2D.h"
#include "Timer.h"

const double pi = 3.1415926535897;

int main() {
  // configures the number of meshes to try (should be between 1 and 8)
  const int num_meshes = 6;
  // configures the number of repetitions computed per mesh to average over
  const int num_tries = 10;

  // source function f
  auto f = [](const Eigen::Vector2d &x) {
    return (8.0 * pi * pi + 1) * std::cos(2 * pi * x(0)) *
           std::cos(2 * pi * x(1));
  };

  // solve the system
  FESolver solver(f);
  FESolver solver_inefficient(f, 1);

  Eigen::Matrix<double, num_meshes, 3> times =
      Eigen::MatrixXd::Zero(num_meshes, 3);

  Timer timer;

  for (int i = 0; i < num_meshes; i++) {
    // load the mesh from file
    TriaMesh2D mesh("meshes/Square" + std::to_string(i + 1) + ".txt");
    // save the number of degrees of freedom
    times(i, 0) = mesh.Coordinates.rows();
    // solve the problem on the mesh num_tries often taking the time for both
    // the efficient and the inefficient method
    Eigen::VectorXd mu;
    for (int j = 0; j < num_tries; j++) {
      std::cout << "Solve " << j + 1 << "/" << num_tries << " of mesh " << i + 1
                << "/" << num_meshes << std::endl;

      // efficient solve
      timer.start();
      mu = solver.Solve(mesh);
      times(i, 1) += timer.elapsed();

      // inefficient solve
      timer.start();
      mu = solver_inefficient.Solve(mesh);
      times(i, 2) += timer.elapsed();
    }
    // average the data and normalize to seconds (measurement is in
    // milliseconds)
    times(i, 1) /= num_tries * 1000.0;
    times(i, 2) /= num_tries * 1000.0;
  }
  // prepare data to be plotted
  mglData x(times.col(0).size(), times.col(0).data());
  mglData y_efficient(times.col(1).size(), times.col(1).data());
  mglData y_inefficient(times.col(2).size(), times.col(2).data());
  // initialize graph and set ranges
  mglGraph graph;
  graph.SetRange('x', x);
  graph.SetRange('y', y_inefficient);
  // configure a log-log plot
  graph.SetFunc("lg(x)", "lg(y)");
  // initialize the axes and set their labels
  graph.SetFontSizePT(8);
  graph.Axis();
  graph.Box();
  graph.Label('x', "Ndof", 0);
  graph.Label('y', "Computational time (in sec.)", 0);
  // plot the data
  graph.Plot(x, y_efficient, "og");
  graph.Plot(x, y_inefficient, "ob");
  // configure the legend
  graph.AddLegend("Efficient matrix assembly", "og");
  graph.AddLegend("Inefficient matrix assembly", "ob");
  graph.Legend(1);
  // write the plot to a file
  graph.WriteEPS("efficiency.eps");
  return 0;
}
