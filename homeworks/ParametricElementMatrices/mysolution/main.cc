
/**
 * @file main.cc
 * @brief NPDE homework ParametricElementMatrices code
 * @author Simon Meierhans
 * @date 14/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "ansiotropic_diffusion_element_matrix_provider.h"
#include "fe_source_elem_vec_provider.h"
#include "impedance_boundary_edge_matrix_provider.h"

int main() {
  // use test mesh to set up finite element space
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1);
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  const lf::base::size_type N_dofs(dofh.NoDofs());

  // interpolate function w(x) = sin(|x|)
  auto w_func = [](Eigen::Vector2d x) -> double { return std::sin(x.norm()); };
  lf::uscalfe::MeshFunctionGlobal mf_w_func{w_func};
  auto w = lf::uscalfe::NodalProjection<double>(fe_space, mf_w_func);

  // initialize Galerkin matrix
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);

  // integrate over domain
  auto d = [](Eigen::Vector2d x) -> Eigen::Vector2d { return x; };
  auto elmat_builder =
      ParametricElementMatrices::AnisotropicDiffusionElementMatrixProvider(d);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // integrate over the boundary
  auto edgemat_builder =
      ParametricElementMatrices::ImpedanceBoundaryEdgeMatrixProvider(fe_space,
                                                                     w);
  lf::assemble::AssembleMatrixLocally(1, dofh, dofh, edgemat_builder, A);

  // initialize load vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();

  // integrate over domain to compute the load vector
  auto elvec_builder =
      ParametricElementMatrices::FESourceElemVecProvider(fe_space, w);
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // solve linear system to obtain the solution
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  Eigen::VectorXd sol_vec = solver.solve(phi);

  std::cout << "Solved system: \n" << sol_vec << "\n";
}
