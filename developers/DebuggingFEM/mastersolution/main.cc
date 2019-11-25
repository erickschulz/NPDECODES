/**
 * @file main.cc
 * @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans
 * @date 27/03/2019
 * @copyright Developed at ETH Zurich
 */

#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <memory>

#include "local_laplace_qfe.h"
#include "qfe_interpolator.h"
#include "qfe_provider_tester.h"

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/refinement/refinement.h>

using size_type = lf::base::size_type;

int main() {
  // read mesh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), CURRENT_SOURCE_DIR "/meshes/square_64.msh");
  auto mesh = reader.mesh();

  // refine mesh
  const int reflevels = 4;
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh, reflevels);
  lf::refinement::MeshHierarchy &multi_mesh{*multi_mesh_p};
  lf::base::size_type L = multi_mesh.NumLevels();

  // vector holding pointers to the different element matrix providers
  std::vector<std::unique_ptr<DebuggingFEM::EntityMatrixProvider>> element_matrix_provider;
  element_matrix_provider.emplace_back(std::make_unique<DebuggingFEM::LocalLaplaceQFE1>());
  element_matrix_provider.emplace_back(std::make_unique<DebuggingFEM::LocalLaplaceQFE2>());
  element_matrix_provider.emplace_back(std::make_unique<DebuggingFEM::LocalLaplaceQFE3>());
  int num_emp = element_matrix_provider.size();

  // vectors accumulating the computed errors
  Eigen::MatrixXd H1SMerr(L, num_emp);
  Eigen::VectorXi N(L);

  for (int level = 0; level < L; ++level) {
    // set up fespace and dof handler for the mesh at the current level
    auto mesh_p = multi_mesh.getMesh(level);
    const lf::assemble::UniformFEDofHandler dofh(mesh_p,
                                           {{lf::base::RefEl::kPoint(), 1},
                                            {lf::base::RefEl::kSegment(), 1},
                                            {lf::base::RefEl::kTria(), 0},
                                            {lf::base::RefEl::kQuad(), 1}});

    // Dimension of finite element space
    N[level] = dofh.NumDofs();

    // compute error for each element matrix provider
    for (int i = 0; i < num_emp; ++i) {
      // function to interpolateOntoQuadFE
      auto f = [](const Eigen::VectorXd &x) -> double {
        return std::exp(x(1) * x(1) + x(0) * x(0));
      };
      // Matrix in triplet format holding Galerkin for LocalLaplaceQFEX matrix,
      // zero initially.
      DebuggingFEM::QFEProviderTester qfe_provider_tester(dofh, *element_matrix_provider[i]);

      // compute the energy / error
#if SOLUTION
      double energy = qfe_provider_tester.energyOfInterpolant(f);
      H1SMerr(level, i) = std::abs(23.76088 - energy);
#else
      //====================
      // Your code goes here
      //====================
#endif

    }

  }

  // Tabular output of the results
  std::cout << std::left << std::setw(10) << "N" << std::setw(20)
            << "Assember 1" << std::setw(20) << "Assembler 2" << std::setw(20)
            << "Assembler 3" << std::endl;
  for (int level = 0; level < L; ++level) {
    std::cout << std::left << std::setw(10) << N[level] << std::setw(20)
              << H1SMerr(level, 0) << std::setw(20) << H1SMerr(level, 1) << std::setw(20)
              << H1SMerr(level, 2) << std::endl;
  }

  /*// Write .csv file and plot it by a python script
  Eigen::MatrixXd data(L, num_emp + 1);
  data.col(0) = N.cast<double>();
  data.rightCols(num_emp) = H1SMerr;
  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream error_file;
  error_file.open("error.csv");
  error_file << data.format(CSVFormat) << std::endl;
  error_file.close();
  std::cout << "Generated " CURRENT_BINARY_DIR "/error.csv" << std::endl;
  std::system("python3 " CURRENT_SOURCE_DIR "/mastersolution/plot_error.py " CURRENT_BINARY_DIR "/error.csv " CURRENT_BINARY_DIR "/error.png");*/

  return 0;
}
