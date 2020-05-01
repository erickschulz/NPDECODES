/**
 * @file nonlinschroedingerequation_main.cc
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 22.04.2020
 * @copyright Developed at ETH Zurich
 */

#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <utility>

#include <Eigen/Core>
#include <Eigen/SparseLU>

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

#include "nonlinschroedingerequation.h"

using namespace std::literals::complex_literals;

int main() {

  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), CURRENT_SOURCE_DIR "/../meshes/square_64.msh");
  auto mesh_p = reader.mesh();

  //lf::io::VtkWriter vtk_writer(mesh_p, "filename.vtk");

  auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<std::complex<double>>>(mesh_p);
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  lf::assemble::COOMatrix<double> D_COO(N_dofs, N_dofs);
  NonLinSchroedingerEquation::MassElementMatrixProvider mass_emp;
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, mass_emp, D_COO);
  Eigen::SparseMatrix<double> D = D_COO.makeSparse();
  Eigen::SparseMatrix<std::complex<double>> Minv = 1i * D;
  Eigen::VectorXcd m = Minv.diagonal();
  Minv.diagonal() = m.cwiseInverse();
  //std::cout << Eigen::MatrixXd(M) << std::endl;


  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  // NonLinSchroedingerEquation::StiffnessElementMatrixProvider
  lf::uscalfe::LinearFELaplaceElementMatrix stiffness_emp;
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, stiffness_emp, A_COO);
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  //std::cout << Eigen::MatrixXd(A) << std::endl;

  int M = 100;
  double T = 1.0;
  double tau = T / M;

  Eigen::SparseMatrix<std::complex<double>> I(N_dofs, N_dofs);
  I.setIdentity();
  Eigen::SparseMatrix<std::complex<double>> MinvA = Minv * A;
  Eigen::SparseMatrix<std::complex<double>> B_plus = I + 0.25 * tau * MinvA;
  Eigen::SparseMatrix<std::complex<double>> B_minus = I - 0.25 * tau * MinvA;

  auto f = [tau] (std::complex<double> z) { return std::exp(-1i * tau * std::norm(z)) * z; };

  Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;
  solver.compute(B_minus);

  //auto u0 = [] (Eigen::Vector2d x) -> std::complex<double> { return std::exp(1i * 3.14159265359 * x.sum()); };  
  //auto u0 = [] (Eigen::Vector2d x) -> std::complex<double> { return -4.0 * std::exp(4.0 * x.sum()) / (std::exp(2.0 * x.sum()) + 1.0); };
  auto u0 = [] (Eigen::Vector2d x) -> std::complex<double> { return std::exp(-x.squaredNorm()); };
  lf::mesh::utils::MeshFunctionGlobal mf_u0{u0};
  Eigen::VectorXcd mu = lf::uscalfe::NodalProjection(*fe_space, mf_u0);
  //std::cout << mu.transpose().format(CSVFormat) << std::endl;
  Eigen::VectorXd norm(M + 1);
  Eigen::VectorXd E_kin(M + 1);
  Eigen::VectorXd E_int(M + 1);
  for (int i = 0; i < M; ++i) {
    // Compute norm and energy along the solution
    norm(i) = NonLinSchroedingerEquation::Norm(mu, D);
    E_kin(i) = NonLinSchroedingerEquation::KineticEnergy(mu, A);
    E_int(i) = NonLinSchroedingerEquation::InteractionEnergy(mu, D);
    // Timestep tau according to Strang splitting
    mu = solver.solve(B_plus * mu);
    mu = mu.unaryExpr(f);
    mu = solver.solve(B_plus * mu);
  }
  norm(M) = NonLinSchroedingerEquation::Norm(mu, D);
  E_kin(M) = NonLinSchroedingerEquation::KineticEnergy(mu, A);
  E_int(M) = NonLinSchroedingerEquation::InteractionEnergy(mu, D);
  //std::cout << mu.transpose().format(CSVFormat) << std::endl;
  //std::cout << norm.transpose().format(CSVFormat) << std::endl;
  //std::cout << (E_kin + E_int).transpose().format(CSVFormat) << std::endl;

  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(M + 1, 0.0, T);

  std::ofstream energies_csv;
  energies_csv.open("energies.csv");
  energies_csv << t.transpose().format(CSVFormat) << std::endl;
  energies_csv << E_kin.transpose().format(CSVFormat) << std::endl;
  energies_csv << E_int.transpose().format(CSVFormat) << std::endl;
  energies_csv.close();

  std::cout << "Generated " CURRENT_BINARY_DIR "/energies.csv" << std::endl;
  std::system("python3 " CURRENT_SOURCE_DIR "/plot_energies.py " CURRENT_BINARY_DIR "/energies.csv " CURRENT_BINARY_DIR "/energies.png");

  std::ofstream norm_csv;
  norm_csv.open("norm.csv");
  norm_csv << t.transpose().format(CSVFormat) << std::endl;
  norm_csv << norm.transpose().format(CSVFormat) << std::endl;
  norm_csv.close();

  std::cout << "Generated " CURRENT_BINARY_DIR "/norm.csv" << std::endl;
  std::system("python3 " CURRENT_SOURCE_DIR "/plot_norm.py " CURRENT_BINARY_DIR "/norm.csv " CURRENT_BINARY_DIR "/norm.png");

  return 0;
}

  /*Eigen::VectorXd quadrature_weights(N_dofs);
  NonLinSchroedingerEquation::TrapezoidalElementVectorProvider trapezoidal_evp;
  lf::assemble::AssembleVectorLocally(0, dofh, trapezoidal_evp, quadrature_weights);
  //std::cout << quadrature_weights.transpose().format(CSVFormat) << std::endl;*/
