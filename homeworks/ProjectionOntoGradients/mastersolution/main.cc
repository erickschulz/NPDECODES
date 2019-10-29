#include <iostream>
#include <memory>

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include "gradprojection.h"

int main() {
  // for this exercise a main file is not required
  // but feel free to use it to call some of your functions for debugging
  // purposes
  // BEGIN_SOLUTION
  const auto f = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x(1), x(0));
  };
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  // Compute solution
  const Eigen::VectorXd sol_vec = ProjectionOntoGradients::projectOntoGradients(dofh, f);
  std::cout << sol_vec << std::endl;

  // END_SOLUTION

  std::cout << "You may use this main file to call your function "
               "ProjectionOntoGradients::ProjectionOntoGradients"
            << std::endl;
}
