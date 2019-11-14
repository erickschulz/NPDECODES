/**
 * @file
 * @brief NPDE homework CoupledSecondOrderBVP
 * @author Erick Schulz
 * @date 13/11/2019
 * @copyright Developed at ETH Zurich
 */

#include "coupledsecondorderbvp.h"

using namespace CoupledSecondOrderBVP;

int main(int /*argc*/, const char** /*argv*/) {
  // Load mesh into a Lehrfem++ object
  std::string mesh_file =
      CURRENT_SOURCE_DIR "/meshes/hex" + std::to_string(0) + ".msh";
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>

  // Load finite element space
  // We discretization by means of piecewise QUADRATIC lagrangian FE
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Solve the coupled boundary value problem
  Eigen::VectorXd sol_vec = solveCoupledBVP(fe_space, 3.0);
}
