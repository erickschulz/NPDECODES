/**
 * @file
 * @brief NPDE homework CoupledSecondOrderBVP
 * @author Amélie Loher
 * @date 01/12/2019
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>

#include "../coupledsecondorderbvp.h"

namespace CoupledSecondOrderBVP::test {

TEST(CoupledSecondOrderBVP, dropMatrixRowsAndColumns) {
    
  std::string mesh_file = 
  		CURRENT_SOURCE_DIR "/meshes/simple.msh";
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>

  // Load finite element space
  // We discretization by means of piecewise QUADRATIC lagrangian FE
  auto fe_space = std::make_shared<CoupledSecondOrderBVP::FeSpaceLagrangeO2<double>>(mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  
  // Obtain an array of boolean flags for the nodes of the mesh, 'true'
  // indicates that the node lies on the boundary
  auto nodes_bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  // Similarly for edges
  auto edges_bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Index predicate for the selectvals FUNCTOR of dropMatrixRowsAndColumns
  auto bd_selector = [&nodes_bd_flags, &edges_bd_flags,
                      &dofh](unsigned int idx) -> bool {
    if (dofh.Entity(idx).RefEl() == lf::base::RefElType::kPoint) {
      return nodes_bd_flags(dofh.Entity(idx));

    } else {
      return edges_bd_flags(dofh.Entity(idx));

    }
  };
  // Coefficients 
  auto const_one = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  auto const_zero = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  lf::assemble::COOMatrix<double> A0(N_dofs, N_dofs);
 
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(const_one), decltype(const_zero)>
      A0_builder(fe_space, const_one, const_zero);

  lf::assemble::AssembleMatrixLocally(0,dofh,dofh,A0_builder,A0);
  // Test dropMatrixRowsAndColumns
  dropMatrixRowsAndColumns(bd_selector, A0);

  const std::vector<Eigen::Triplet<double>> A0_triplets_vec = A0.triplets();
  // Matrix Storing all nnz entries of A0  
  Eigen::MatrixXd A0triplets;

  for (auto &triplet : A0_triplets_vec) {
    A0triplets << triplet.value(); 
  }
  // Reference values for nnz entries of A0
  Eigen::MatrixXd refval;
  refval << 1, -0.666667, -0.666667, -0.666667, 2.66667, -6.10623e-16, -0.666667, -6.10623e-16,
				2.66667, 1, -0.666667, -0.666667, -0.666667, 2.66667, -6.10623e-16, -0.666667,
				-6.10623e-16, 2.66667, 1, -0.666667, -0.666667,-0.666667,2.66667,-6.10623e-16,
				-0.666667,-6.10623e-16,2.66667,1,-0.666667,-0.666667,-0.666667,2.66667,
				-6.10623e-16,-0.666667,-6.10623e-16,2.66667,1,1,1,
				1,1,1,1,1;

  double tol = 1.0e-8;
  ASSERT_NEAR(refval.norm(), A0triplets.norm(), tol);

}


TEST(CoupledSecondOrderBVP, dropMatrixRows) {

}


TEST(CoupledSecondOrderBVP, solveCoupledBVP) {


}

} // end namespace
