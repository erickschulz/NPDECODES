/**
 * @file nonlinschroedingerequation_main.cc
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 22.04.2020
 * @copyright Developed at ETH Zurich
 */

#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <utility>

#include <Eigen/Core>
#include <Eigen/SparseCholesky>

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

#include "nonlinschroedingerequation.h"

using namespace std::literals::complex_literals;

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
  auto inv = [] (double x) { return x == 0.0 ? 0.0 : 1.0 / x; };
  Eigen::SparseMatrix<std::complex<double>> Minv = -1i * D.unaryExpr(inv);
  //std::cout << Eigen::MatrixXcd(Minv) << std::endl;

  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  NonLinSchroedingerEquation::StiffnessElementMatrixProvider stiffness_emp;
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, stiffness_emp, A_COO);
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  //std::cout << Eigen::MatrixXd(A) << std::endl;

  int M = 100;
  double T = 1.0;
  double tau = T / M;

  Eigen::SparseMatrix<std::complex<double>> I(N_dofs, N_dofs);
  I.setIdentity();
  Eigen::SparseMatrix<std::complex<double>> temp = tau * Minv * A.cast<std::complex<double>>();
  Eigen::SparseMatrix<std::complex<double>> B_plus = I + temp;
  Eigen::SparseMatrix<std::complex<double>> B_minus = I - temp;
/*
  auto f = [tau] (std::complex<double> z) { return -1i * z * std::norm(z) * tau; }

  SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver;
  Eigen::VectorXcd mu;
  for (int i = 0; i < M; ++i) {
    mu = mu.unaryExpr(f)
  }
*/
  return 0;
}
