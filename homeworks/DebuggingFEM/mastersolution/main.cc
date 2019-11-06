/**
 * @file main.cc
 * @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans
 * @date 27/03/2019
 * @copyright Developed at ETH Zurich
 */

#include <functional>
#include <iomanip>

#include "local_laplace_qfe.h"
#include "qfe_interpolator.h"
#include "qfe_provider_tester.h"

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/refinement/refinement.h>

#include <boost/filesystem.hpp>

using size_type = lf::base::size_type;

int main()
{
  // read mesh
  boost::filesystem::path here = __FILE__;
  auto square_path = here.parent_path().parent_path() / "meshes/square_64.msh";

  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), square_path.string());
  auto mesh = reader.mesh();

  // refine mesh
  const int reflevels = 4;
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh, reflevels);
  lf::refinement::MeshHierarchy &multi_mesh{*multi_mesh_p};
  size_type L = multi_mesh.NumLevels();

  // function to interpolateOntoQuadFE
  auto f = [](Eigen::VectorXd x) -> double {
    return std::exp(x(1) * x(1) + x(0) * x(0));
  };

  // vectors accumulating the computed errors
  std::vector<double> H1SMerr_1(L);
  std::vector<double> H1SMerr_2(L);
  std::vector<double> H1SMerr_3(L);
  std::vector<double> ndofs;

  for (int level = 0; level < L; ++level)
  {
    // set up fespace and dof handler for the mesh at the current level
    auto mesh_p = multi_mesh.getMesh(level);
    lf::assemble::UniformFEDofHandler dofh(mesh_p,
                                           {{lf::base::RefEl::kPoint(), 1},
                                            {lf::base::RefEl::kSegment(), 1},
                                            {lf::base::RefEl::kTria(), 0},
                                            {lf::base::RefEl::kQuad(), 1}});
    // Dimension of finite element space`
    const size_type N_dofs(dofh.NumDofs());
    // Matrix in triplet format holding Galerkin for LocalLaplaceQFE1 matrix,
    // zero initially.
    lf::assemble::COOMatrix<double> A_1(N_dofs, N_dofs);
    auto element_matrix_provider_1 = DebuggingFEM::LocalLaplaceQFE1();
    DebuggingFEM::QFEProviderTester qfe_provider_tester_1(
        dofh, element_matrix_provider_1);

    // Matrix in triplet format holding Galerkin for LocalLaplaceQFE2 matrix,
    // zero initially.
    lf::assemble::COOMatrix<double> A_2(N_dofs, N_dofs);
    auto element_matrix_provider_2 = DebuggingFEM::LocalLaplaceQFE2();
    DebuggingFEM::QFEProviderTester qfe_provider_tester_2(
        dofh, element_matrix_provider_2);

    // Matrix in triplet format holding Galerkin for LocalLaplaceQFE3 matrix,
    // zero initially.
    lf::assemble::COOMatrix<double> A_3(N_dofs, N_dofs);
    auto element_matrix_provider_3 = DebuggingFEM::LocalLaplaceQFE3();
    DebuggingFEM::QFEProviderTester qfe_provider_tester_3(
        dofh, element_matrix_provider_3);

    // interpolateOntoQuadFE exact solution
    Eigen::VectorXd u_exact = DebuggingFEM::interpolateOntoQuadFE(dofh, f);

    ndofs.push_back(N_dofs);

    // compute the energies
    /* SOLUTION_BEGIN */
    double energy_1 = qfe_provider_tester_1.energyOfInterpolant(f);
    double h1s_error_1 = std::abs(23.76088 - energy_1);
    H1SMerr_1[level] = h1s_error_1;

    double energy_2 = qfe_provider_tester_2.energyOfInterpolant(f);
    double h1s_error_2 = std::abs(23.76088 - energy_2);
    H1SMerr_2[level] = h1s_error_2;

    double energy_3 = qfe_provider_tester_3.energyOfInterpolant(f);
    double h1s_error_3 = std::abs(23.76088 - energy_3);
    H1SMerr_3[level] = h1s_error_3;
    /* SOLUTION_END */
  }

  // Tabular output of the results
  std::cout << std::left << std::setw(10) << "N" << std::setw(20)
            << "Assember 1" << std::setw(20) << "Assembler 2" << std::setw(20)
            << "Assembler 3" << std::endl;
  for (int l = 0; l < L; ++l)
  {
    std::cout << std::left << std::setw(10) << ndofs[l] << std::setw(20)
              << H1SMerr_1[l] << std::setw(20) << H1SMerr_2[l] << std::setw(20)
              << H1SMerr_3[l] << std::endl;
  }

  return 0;
}
