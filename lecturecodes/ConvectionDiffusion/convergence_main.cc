/**
 * @file upwindquadrature_main.cc
 * @brief Convergence study for 3 methods to solve CD BVP
 * @author Philippe Peter
 * @date July 2020
 * @copyright Developed at SAM, ETH Zurich
 */
#include <lf/fe/fe.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "cd_tools.h"
#include "standard_fem.h"
#include "supg.h"
#include "upwind.h"

int main() {
  // parameter functions:
  // velocity field
  const Eigen::Vector2d v(2.0, 3.0);
  const auto velocity = [&v](const Eigen::Vector2d& x) { return v; };
  // diffusion coefficient
  const double eps = 1.0;
  const auto epsilon = [&eps](const Eigen::Vector2d& x) { return eps; };
  // exact solution:
  const auto u_exact = [&eps](const Eigen::Vector2d& x) {
    return x(0) * x(1) * x(1) - x(1) * x(1) * std::exp(2 * (x(0) - 1) / eps) -
           x(0) * std::exp(3 * (x(1) - 1) / eps) +
           std::exp(2 * (x(0) - 1) / eps + 3 * (x(1) - 1) / eps);
  };
  // gradient of the exact solution
  const auto u_grad_exact = [&eps](const Eigen::Vector2d& x) {
    Eigen::Vector2d res;
    res(0) = x(1) * x(1) -
             x(1) * x(1) * 2.0 / eps * std::exp(2 * (x(0) - 1) / eps) -
             std::exp(3 * (x(1) - 1) / eps) +
             2.0 / eps * std::exp(2 * (x(0) - 1) / eps + 3 * (x(1) - 1) / eps);
    res(1) = 2 * x(0) * x(1) - 2 * x(1) * std::exp(2 * (x(0) - 1) / eps) -
             x(0) * 3.0 / eps * std::exp(3 * (x(1) - 1) / eps) +
             3.0 / eps * std::exp(2 * (x(0) - 1) / eps + 3 * (x(1) - 1) / eps);
    return res;
  };
  // Laplacian of the exact solution
  const auto u_laplace_exact = [&eps](const Eigen::Vector2d& x) {
    double uxx =
        -4.0 / (eps * eps) * x(1) * x(1) * std::exp(2 * (x(0) - 1) / eps) +
        4.0 / (eps * eps) *
            std::exp(2 * (x(0) - 1) / eps + 3 * (x(1) - 1) / eps);
    double uyy = 2.0 * x(0) - 2.0 * std::exp(2 * (x(0) - 1) / eps) -
                 x(0) * 9.0 / (eps * eps) * std::exp(3 * (x(1) - 1) / eps) +
                 9.0 / (eps * eps) *
                     std::exp(2 * (x(0) - 1) / eps + 3 * (x(1) - 1) / eps);
    return uxx + uyy;
  };
  // source function
  const auto f = [&eps, &v, &u_grad_exact,
                  &u_laplace_exact](const Eigen::Vector2d& x) {
    return -eps * u_laplace_exact(x) + v.transpose() * u_grad_exact(x);
  };

  // boundary conditions
  const auto g = [&u_exact](const Eigen::Vector2d& x) { return u_exact(x); };

  // Construct mesh hierarchy:
  std::unique_ptr<lf::mesh::MeshFactory> top_mesh_factory_ptr =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::utils::TPTriagMeshBuilder builder(std::move(top_mesh_factory_ptr));
  builder.setBottomLeftCorner(Eigen::Vector2d{0.0, 0.0})
      .setTopRightCorner(Eigen::Vector2d{1.0, 1.0})
      .setNumXCells(2)
      .setNumYCells(2);
  std::shared_ptr<lf::mesh::Mesh> top_mesh = builder.Build();

  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(top_mesh, 6);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  multi_mesh.PrintInfo(std::cout);

  // get number of levels:
  unsigned L = multi_mesh.NumLevels();

  // Output file
  std::string file_name = "results_errors.txt";
  std::ofstream file;
  file.open(file_name);
  file << "h, $L^2$-Error u (FEM), $L^2$-Error u (Upwind), $L^2$-Error u "
          "(SUPG) \n";

  // Perform computations on all levels and compute  L2-errors
  for (unsigned l = 0; l < L; ++l) {
    // extract mesh and construct FE space.
    std::shared_ptr<const lf::mesh::Mesh> mesh_p{multi_mesh.getMesh(l)};
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

    // compute solutions using Standard FE, Upwind, SUPG method
    Eigen::VectorXd sol_standard = ConvectionDiffusion::SolveCDBVPStandardFem(
        fe_space, epsilon, velocity, f, g);
    lf::fe::MeshFunctionFE sol_standard_mf(fe_space, sol_standard);

    Eigen::VectorXd sol_stable = ConvectionDiffusion::SolveCDBVPUpwind(
        fe_space, epsilon, velocity, f, g);
    lf::fe::MeshFunctionFE sol_upwind_mf(fe_space, sol_stable);

    Eigen::VectorXd sol_supg =
        ConvectionDiffusion::SolveCDBVPSupg(fe_space, epsilon, velocity, f, g);
    lf::fe::MeshFunctionFE sol_supg_mf(fe_space, sol_supg);

    // Wrap exact solution into mesh function for error computations
    auto mf_solution = lf::mesh::utils::MeshFunctionGlobal(u_exact);

    // Calculate L2 errors:
    double L2err_standard = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::mesh::utils::squaredNorm(sol_standard_mf - mf_solution),
        10));
    double L2err_upwind = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::mesh::utils::squaredNorm(sol_upwind_mf - mf_solution),
        10));
    double L2err_supg = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::mesh::utils::squaredNorm(sol_supg_mf - mf_solution), 10));

    // output
    std::cout << "Level " << l
              << "(h= " << ConvectionDiffusion::MeshWidth(mesh_p)
              << "): " << L2err_standard << ", " << L2err_upwind << ", "
              << L2err_supg << std::endl;
    file << ConvectionDiffusion::MeshWidth(mesh_p) << ", " << L2err_standard
         << ", " << L2err_upwind << ", " << L2err_supg << std::endl;
  }
  file.close();

  // Plot
  std::system("python3 " CURRENT_SOURCE_DIR
              "/plot_convergence.py " CURRENT_BINARY_DIR);

  return 0;
}