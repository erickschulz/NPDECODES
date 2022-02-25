/**
 * @ file
 * @ brief NPDE homework on Nitsche's method
 * @ author R. Hiptmair
 * @ date July 2021
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <cmath>
#include <iomanip>

#include "nitschemethod.h"

int main(int /*argc*/, char** /*argv*/) {
  std::cout << "Homework problem on Nitsche's method" << std::endl;

  // Exact solution u, a harmonic function, will also supply Dirichlet data
  auto u = [](Eigen::Vector2d x) -> double {
    return x[0] * x[0] - x[1] * x[1];
  };
  // Has to be wrapped into a mesh function for error computation
  lf::mesh::utils::MeshFunctionGlobal mf_u{u};

  // Gradient of exact solution
  auto grad_u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return ((Eigen::Vector2d() << 2.0 * x[0], -2.0 * x[1]).finished());
  };
  // Convert into mesh function to use for error computation
  lf::mesh::utils::MeshFunctionGlobal mf_grad_u{grad_u};

  // Obtain a triangular mesh of the unit square from the collection of
  // test meshes
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);

  // Generate a sequence of meshes by regular refinement.
  const int reflevels = 5;
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p,
                                                              reflevels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  // Number of levels
  int L = multi_mesh.NumLevels();

  // Vector for keeping error norms
  std::vector<std::tuple<int, double, double>> errs{};

  // Penalty parameter on coarsest mesh
  double c = 60.0;

  // LEVEL LOOP: Do computations on all levels
  for (int level = 0; level < L; ++level, c *= 3) {
    // Pointer to current mesh
    mesh_p = multi_mesh.getMesh(level);
    // Pointer to lowest-degree Lagrangian finite element space
    auto fes_p =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Assemble Galerkin matrix
    Eigen::SparseMatrix<double> A =
        NitscheMethod::computeNitscheGalerkinMatrix(fes_p, c);
    // Assemble right-hand side vector
    Eigen::VectorXd phi = NitscheMethod::computeNitscheLoadVector(fes_p, u, c);
    // Solve linear system using Eigen's sparse direct elimination
    // Examine return status of solver in case the matrix is singular
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    Eigen::VectorXd mu = solver.solve(phi);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

    // Compute error norms
    // create *mesh_p functions representing solution / gradient of solution
    const lf::fe::MeshFunctionFE mf_sol(fes_p, mu);
    const lf::fe::MeshFunctionGradFE mf_grad_sol(fes_p, mu);
    // compute errors with 3rd order quadrature rules, which is sufficient for
    // piecewise linear finite elements
    double L2err =  // NOLINT
        std::sqrt(lf::fe::IntegrateMeshFunction(
            *mesh_p, lf::mesh::utils::squaredNorm(mf_sol - mf_u), 2));
    double H1serr = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        *mesh_p, lf::mesh::utils::squaredNorm(mf_grad_sol - mf_grad_u), 2));

    std::cout << "Mesh (" << mesh_p->NumEntities(0) << " cells, "
              << mesh_p->NumEntities(1) << " edges, " << mesh_p->NumEntities(2)
              << " nodes): L2err = " << L2err << ", H1serr = " << H1serr
              << std::endl;
    errs.emplace_back(mesh_p->NumEntities(0), L2err, H1serr);
  }

  std::cout << std::left << std::setw(16) << "#cells" << std::right
            << std::setw(16) << "L2 error" << std::setw(16) << "H1 error"
            << std::endl;
  for (const auto& err : errs) {
    auto [N, l2err, h1serr] = err;
    std::cout << std::left << std::setw(16) << N << std::left << std::setw(16)
              << l2err << std::setw(16) << h1serr << std::endl;
  }

  return 0;
}
