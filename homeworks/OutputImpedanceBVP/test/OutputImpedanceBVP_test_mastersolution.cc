/**
 * @file
 * @brief UNITESTS for NPDE homework ErrorEstimatesForTraces
 * @author Erick Schulz
 * @date 28/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "../mastersolution/OutputImpedanceBVP.h"
#include "../mastersolution/evalclass.h"
#include <gtest/gtest.h>

namespace OutputImpedanceBVP::test
{

TEST(OutputImpedanceBVP, computeApproxSolDirichlet)
{
  std::cout << "TESTS - OutputImpedanceBVP " << std::endl;

  // Load mesh into a Lehrfem++ object
  boost::filesystem::path here = __FILE__;
  std::string filename = "/meshes/unitsquare.msh";
  auto mesh_path = here.parent_path().parent_path() / filename;
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_path.string());
  auto mesh_p = reader.mesh(); // type shared_ptr< const lf::mesh::Mesh>

  // Finite element space
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Exact solution and Dirichlet boundary conditions
  Eigen::Vector2d g;
  g << 1.0, 3.0;
  auto uExact = [&g](coord_t x) -> double { return g.dot(x); };
  auto uExact_vec = interpolateData<std::function<double(Eigen::Vector2d)>>(
      fe_space_p, std::move(uExact));

  // Solve BVP
  Eigen::VectorXd uApprox_vec = solveImpedanceBVP(fe_space_p, g);

  ASSERT_TRUE(uApprox_vec.isApprox(uExact_vec));
}

} // namespace OutputImpedanceBVP::test
