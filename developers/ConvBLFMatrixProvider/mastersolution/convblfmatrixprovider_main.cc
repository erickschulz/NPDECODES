/**
 * @file
 * @brief demonstration of assembly of Galerkin linear system in LehrFEM++
 * @author Ralf Hiptmair
 * @date   April 2021
 * @copyright Developed at ETH Zurich
 */

#include "convblfmatrixprovider.h"

int main(int /*argc*/, char** /*argv*/) {
  // Obtain a purely triangular mesh from the collection of LehrFEM++'s
  // built-in meshes
  std::shared_ptr<lf::mesh::Mesh> mesh_p{
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3)};
  // vectors for testing
  const Eigen::Vector2d a(1.0, 2.0);
  const Eigen::Vector2d b(3.0, 2.0);
  // Run test computation
  double itg = cblfdemo::testCDBLF(mesh_p, a, b);
  std::cout << "Integral = " << itg << std::endl;
  return 0;
}
