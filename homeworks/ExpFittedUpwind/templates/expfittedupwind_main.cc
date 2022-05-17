/**
 * @file expfittedupwind_main.cc
 * @brief NPDE homework ExpFittedUpwind
 * @author Amélie Loher, Philippe Peter
 * @date 07.01.2020
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>

#include "expfittedupwind.h"

int main() {
  // Define Mesh-independent Data:
  //====================
  // Your code goes here
  //====================

  // Output file
  std::ofstream L2output;
  L2output.open(CURRENT_BINARY_DIR "/L2error.txt");
  L2output << "No. of dofs, L2 error" << std::endl;

  // generate a mesh hierarchy:
  unsigned int reflevels = 6;
  std::unique_ptr<lf::mesh::MeshFactory> mesh_factory_ptr =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::utils::TPTriagMeshBuilder builder(std::move(mesh_factory_ptr));
  builder.setBottomLeftCorner(Eigen::Vector2d{0.0, 0.0})
      .setTopRightCorner(Eigen::Vector2d{1.0, 1.0})
      .setNumXCells(2)
      .setNumYCells(2);
  auto top_mesh = builder.Build();

  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(top_mesh,
                                                              reflevels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  multi_mesh.PrintInfo(std::cout);

  // get number of levels:
  auto L = multi_mesh.NumLevels();

  // perform computations on all levels:
  for (int l = 0; l < L; ++l) {
    // Compute finite element solution and compute L2 error on current level:
    double L2_err = 1.0;

    // get current mesh and fe space
    auto mesh_p = multi_mesh.getMesh(l);
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

    const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
    const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

    //====================
    // Your code goes here
    //====================

    L2output << N_dofs << ", " << L2_err << std::endl;
    std::cout << N_dofs << "," << L2_err << std::endl;
  }

  L2output.close();

  // Plot the computed L2 error
  std::system("python3 " CURRENT_SOURCE_DIR "/plot_error.py " CURRENT_BINARY_DIR
              "/L2error.txt " CURRENT_BINARY_DIR "/results.eps");

  return 0;
}
