/**
 * @file nonlinschroedingerequation_main.cc
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 22.04.2020
 * @copyright Developed at ETH Zurich
 */

#include <iostream>
#include <memory>
#include <utility>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

#include "nonlinschroedingerequation.h"

int main() {

  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), CURRENT_SOURCE_DIR "/../meshes/square_64.msh");
  auto mesh_p = reader.mesh();

  auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  lf::assemble::COOMatrix<double> D_COO(N_dofs, N_dofs);
  NonLinSchroedingerEquation::MassElementMatrixProvider mass_emp;
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, mass_emp, D_COO);
  Eigen::SparseMatrix<double> D = D_COO.makeSparse();
  //std::cout << Eigen::MatrixXd(D) << std::endl;

  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  NonLinSchroedingerEquation::StiffnessElementMatrixProvider stiffness_emp;
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, stiffness_emp, A_COO);
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  //std::cout << Eigen::MatrixXd(A) << std::endl;

  return 0;
}
