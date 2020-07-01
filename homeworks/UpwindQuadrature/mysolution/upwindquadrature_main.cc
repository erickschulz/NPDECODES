/**
 * @file upwindquadrature_main.cc
 * @brief NPDE homework template main
 * @author Philippe Peter
 * @date June 2020
 * @copyright Developed at SAM, ETH Zurich
 */
#include <cmath>
#include <memory>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include "upwindquadrature.h"

int main() {
  // PARAMETERS
  // mesh specification (number of cells in both sides of the tensor-product
  // triangular mesh)
  int M = 49;

  // coefficient functions:
  // Dirichlet functor
  const auto g = [](const Eigen::Vector2d &x) {
    return x(1) == 0 ? 0.5 - std::abs(x(0) - 0.5) : 0.0;
  };
  lf::mesh::utils::MeshFunctionGlobal mf_g{g};

  // velocity field
  const auto v = [](const Eigen::Vector2d &x) {
    return (Eigen::Vector2d() << -x(1), x(0)).finished();
  };

  // diffusion coefficient
  const double eps = 1e-4;
  lf::mesh::utils::MeshFunctionConstant mf_eps{eps};

  // MESH CONSTRUCTION
  // construct a triangular tensor product mesh on the unit square
  std::unique_ptr<lf::mesh::MeshFactory> mesh_factory_ptr =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::hybrid2d::TPTriagMeshBuilder builder(std::move(mesh_factory_ptr));
  builder.setBottomLeftCorner(Eigen::Vector2d{0.0, 0.0})
      .setTopRightCorner(Eigen::Vector2d{1.0, 1.0})
      .setNumXCells(M)
      .setNumYCells(M);
  std::shared_ptr<lf::mesh::Mesh> mesh_p = builder.Build();

  // DOF HANDLER & FINITE ELEMENT SPACE
  // Construct dofhanlder for linear finite elements on the mesh.
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  // PREPARING DATA TO IMPOSE DIRICHLET CONDITIONS
  // Obtain specification for shape functions on edges
  std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<double>>
      rsf_edge_p = fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());

  // Create a dataset of boolean flags indicating edges on the boundary of the
  // mesh
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

  // Fetch flags and values for degrees of freedom located on Dirichlet
  // boundary.
  auto ess_bdc_flags_values{lf::uscalfe::InitEssentialConditionFromFunction(
      dofh, *rsf_edge_p, bd_flags, mf_g)};

  //============================================================================
  // SOLVE LAPLACIAN WITH NON-HOMOGENEOUS DIRICHLET BC (STANDARD: UNSTABLE)
  //============================================================================
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(dofh.NumDofs(), dofh.NumDofs());

  // ASSEMBLE GALERKIN MATRIX
  // First the part corresponding to the laplacian
  lf::uscalfe::ReactionDiffusionElementMatrixProvider laplacian_provider(
      fe_space, mf_eps, lf::mesh::utils::MeshFunctionConstant(0.0));
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, laplacian_provider, A);

  // Next part corresponding to the convection term:
  UpwindQuadrature::ConvectionElementMatrixProvider convection_provider(v);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, convection_provider, A);

  // RIGHT-HAND SIDE VECTOR
  Eigen::VectorXd phi(dofh.NumDofs());
  phi.setZero();

  // IMPOSE DIRICHLET CONDITIONS:
  // Eliminate Dirichlet dofs from linear system
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&ess_bdc_flags_values](lf::uscalfe::glb_idx_t gdof_idx) {
        return ess_bdc_flags_values[gdof_idx];
      },
      A, phi);

  // SOLVE LINEAR SYSTEM
  Eigen::SparseMatrix A_crs = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  Eigen::VectorXd sol_vec = solver.solve(phi);

  // OUTPUT RESULTS TO VTK FILe
  // construct mesh function representing finite element solution
  lf::uscalfe::MeshFunctionFE mf_sol(fe_space, sol_vec);
  // construct vtk writer
  lf::io::VtkWriter vtk_writer(
      mesh_p, CURRENT_BINARY_DIR "/upwind_quadrature_solution_unstable.vtk");
  // output data
  vtk_writer.WritePointData("upwind_quadrature_solution_unstable", mf_sol);

  //============================================================================
  // SOLVE LAPLACIAN WITH NON-HOMOGENEOUS DIRICHLET BC (UPWIND: STABLE)
  //============================================================================
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A_stable(dofh.NumDofs(), dofh.NumDofs());

  // ASSEMBLE GALERKIN MATRIX
  // First the part corresponding to the laplacian, computed using standard
  // Galerkin approach
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, laplacian_provider,
                                      A_stable);

  // Next part corresponding to the convection term, computed using upwind
  // quadrature:
  UpwindQuadrature::UpwindConvectionElementMatrixProvider
      convection_provider_stable(v, UpwindQuadrature::initializeMasses(mesh_p));
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, convection_provider_stable,
                                      A_stable);

  // RIGHT-HAND SIDE VECTOR
  Eigen::VectorXd phi_stable(dofh.NumDofs());
  phi_stable.setZero();

  // IMPOSE DIRICHLET CONDITIONS:
  // Eliminate Dirichlet dofs from linear system
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&ess_bdc_flags_values](lf::uscalfe::glb_idx_t gdof_idx) {
        return ess_bdc_flags_values[gdof_idx];
      },
      A_stable, phi_stable);

  // SOLVE LINEAR SYSTEM
  Eigen::SparseMatrix A_stable_crs = A_stable.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_stable;
  solver_stable.compute(A_stable_crs);
  Eigen::VectorXd sol_vec_stable = solver_stable.solve(phi_stable);

  // OUTPUT RESULTS TO VTK FILe
  // construct mesh function representing finite element solution
  lf::uscalfe::MeshFunctionFE mf_sol_stable(fe_space, sol_vec_stable);
  // construct vtk writer
  lf::io::VtkWriter vtk_writer_stable(
      mesh_p, CURRENT_BINARY_DIR "/upwind_quadrature_solution_stable.vtk");
  // output data
  vtk_writer_stable.WritePointData("upwind_quadrature_solution_stable",
                                   mf_sol_stable);

  return 0;
}
