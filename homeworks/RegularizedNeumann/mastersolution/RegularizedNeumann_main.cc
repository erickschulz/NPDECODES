/**
 * @ file RegularizedNeumann_main.cc
 * @ brief NPDE homework RegularizedNeumann code
 * @ author Christian Mitsch
 * @ date 11.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "regNeumann.h"

int main() {
  /* BEGIN_SOLUTION */
  const auto f = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  const auto h = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 3.0);

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Compute solution
  auto result_c = RegularizedNeumann::getGalerkinLSE_dropDof(fe_space, f, h);
  auto result_f = RegularizedNeumann::getGalerkinLSE_augment(fe_space, f, h);
  /* END_SOLUTION */
  std::cout << "You can use the mainfile to call your functions" << std::endl;
}
