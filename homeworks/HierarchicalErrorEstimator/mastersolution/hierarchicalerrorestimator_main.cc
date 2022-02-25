/**
 * @ file
 * @ brief NPDE homework on hierarchical error estimation
 * @ author Ralf Hiptmair
 * @ date July 2021
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <lf/mesh/test_utils/test_meshes.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "hierarchicalerrorestimator.h"

int main(int /*argc*/, char** /*argv*/) {
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
  std::vector<std::tuple<int, double, double, double>> errs{};

  // LEVEL LOOP: Do computations on all levels
  for (int level = 0; level < L; ++level) {
    mesh_p = multi_mesh.getMesh(level);
    auto [L2err, H1serr, ree] = HEST::solveAndEstimate(mesh_p);
    errs.emplace_back(mesh_p->NumEntities(0), L2err, H1serr, ree);
  }

  std::cout << std::left << std::setw(16) << "#cells" << std::right
            << std::setw(16) << "L2 error" << std::setw(16) << "H1 error"
            << std::setw(16) << "Estimate" << std::endl;
  for (const auto& err : errs) {
    auto [N, l2err, h1serr, est] = err;
    std::cout << std::left << std::setw(16) << N << std::left << std::setw(16)
              << l2err << std::setw(16) << h1serr << std::setw(16) << est
              << std::endl;
  }
  return 0;
}
