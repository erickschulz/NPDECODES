/**
 * @file regularizedneumann_main.cc
 * @brief NPDE homework RegularizedNeumann code
 * @author Christian Mitsch, Philippe PEter
 * @date March 2020
 * @copyright Developed at ETH Zurich
 */

#include <iostream>
#include <memory>

#include <Eigen/Core>

#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include "regularizedneumann.h"

int main() {

  std::cout << "You can use the mainfile to call your functions" << std::endl;

  const auto f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  const auto h = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 3.0);

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Compute solution
  auto result_c = RegularizedNeumann::getGalerkinLSE_dropDof(fe_space, f, h);
  auto result_f = RegularizedNeumann::getGalerkinLSE_augment(fe_space, f, h);

  return 0;
}
