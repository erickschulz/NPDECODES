#include "test_problem.h"

namespace TestProblem {

// A trivial function
int add(int a, int b) {
  int r;

  /* SOLUTION_BEGIN */
  r = a + b;
  /* SOLUTION_END */

  return r;
}

// LehrFEM++ function calculating the shape regularity measure
double shapeRegularityMeasure(const lf::mesh::Mesh &mesh) {
  // Dimension of the meshed manifold
  const lf::base::size_type dim = mesh.DimMesh();
  // shape regularity measure
  double rho_mesh = 0.0;
  // Run through all the cells of the mesh
  for (const lf::mesh::Entity &cell : mesh.Entities(0)) {
    // Topological type of cell
    const lf::base::RefEl ref_el{cell.RefEl()};
    // Obtain geometry of the current cell
    const lf::geometry::Geometry *geo_p{cell.Geometry()};
    LF_ASSERT_MSG(geo_p != nullptr, "Invalid geometry for " << cell);
    // Area of current cell
    const double area = lf::geometry::Volume(*geo_p);
    // Obtain corner coordinates
    Eigen::MatrixXd vertices{lf::geometry::Corners(*geo_p)};
    // Number of vertices
    const lf::base::size_type num_nodes{ref_el.NumNodes()};
    LF_ASSERT_MSG(vertices.cols() == num_nodes, "num nodes mismatch");
    // Compute length of longest edge
    double hK = 0.0;
    for (int l = 0; l < num_nodes; ++l) {
      auto edge_vec{vertices.col((l + 1) % num_nodes) - vertices.col(l)};
      hK = std::max(hK, edge_vec.norm());
    }
    // Local shape regularity measure
    const double rho_K = (std::pow(hK, dim) / area);
    // Global shape regularity measure
    rho_mesh = std::max(rho_mesh, rho_K);
  }
  return rho_mesh;
}

}  // namespace TestProblem
