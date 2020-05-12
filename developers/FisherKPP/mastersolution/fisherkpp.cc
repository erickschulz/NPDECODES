/** @file fisherkpp.cc
 *  @brief Homework Problem Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 12.05.20
 *  @copyright Developed at SAM, ETH Zurich
 */

#include <memory>
#include <utility>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

namespace FisherKPP {


/* Function for the assembly of both Galerkin Matrices, 
 * the Mass matrix and the Stiffness matrix.
 */
template <typename DIFF_COEFF>
std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
  assembleGalerkinMatrices(const lf::assemble::DofHandler& dofh, DIFF_COEFF &&c) {

  std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> A_M;
  
  // Obtain mesh and finite element space
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofh.Mesh();
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  
  // Matrix format for assembly 
  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  lf::assemble::COOMatrix<double> M_COO(N_dofs, N_dofs);
  
  // Function handles for Stiffness and Mass matrix
  auto one = [] (Eigen::Vector2d x) -> double { return 1.0; };
  auto zero = [] (Eigen::Vector2d x) -> double { return 0.0; };
  
  // Mesh functions
  lf::mesh::utils::MeshFunctionGlobal mf_c{c};
  lf::mesh::utils::MeshFunctionGlobal mf_one{one};
  lf::mesh::utils::MeshFunctionGlobal mf_zero{zero};
  
  // Element matrix provider
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(mf_c), decltype(mf_zero)> elMat_Stiff(fe_space, mf_c, mf_zero);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(mf_zero), decltype(mf_one)> elMat_Mass(fe_space, mf_zero, mf_one);
  
  // Assembly
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elMat_Stiff, A_COO);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elMat_Mass, M_COO);
  
  // Sparse matrix format
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  Eigen::SparseMatrix<double> M = M_COO.makeSparse();

  A_M = std::make_pair(A, M);
  return A_M;

}


class StrangSplit {

  public:
    // Disabled constructors
    StrangSplit() = delete;
    StrangSplit(const StrangSplit &) = delete;
    StrangSplit(StrangSplit &&) = delete;
    StrangSplit &operator=(const StrangSplit &) = delete;
    StrangSplit &operator=(const StrangSplit &&) = delete;
    // Main constructor 
    template <typename DIFF_COEFF>
    explicit StrangSplit(const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space,
                         double T, unsigned int m, double lambda, DIFF_COEFF &&c);

    // Destructor 
    virtual ~StrangSplit() = default;

	/* Member Function StrangSplit
     * Computes the Evolution Operator for the linear parabolic diffusion term
     */
    Eigen::VectorXd diffusionEvolutionOperator(double tau, const Eigen::VectorXd &mu) {
	  Eigen::VectorXd evol_op;

      // Precomputation 
      solver.compute(M_ + tau * xi_ * A_);
      LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
      Eigen::VectorXd rhs = -A_ * mu;

      // First stage SDIRK-2 
      Eigen::VectorXd k1 = solver.solve(rhs);
      LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
      // Second stage SDIRK-2 
      Eigen::VectorXd k2 = solver.solve(rhs - tau * (1 - xi_) * A_ * k1);
      LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
      // Recover solution
      evol_op = mu + tau * (1 - xi_) * k1 + tau * xi_ * k2;
      return evol_op;
    }
    
    /* Member Function StrangSplit
     * Computes the Evolution for m_ timesteps 
     */
	Eigen::VectorXd Evolution(const Eigen::VectorXd &cap, const Eigen::VectorXd &mu) {
	  // Obtain dofhandler
      const lf::assemble::DofHandler &dofh{fe_space_->LocGlobMap()};
      const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  	  Eigen::VectorXd sol(N_dofs);
  	  Eigen::VectorXd sol_cur(N_dofs);
      Eigen::VectorXd sol_next(N_dofs);

      Eigen::VectorXd ones = Eigen::VectorXd::Ones(N_dofs);

  	  // step size
  	  double tau = T_ / m_;

  	  /* Strang Splitting Method 
   	   * First half time step [0, tau/2]: Diffusion
       */
  	  sol_cur = diffusionEvolutionOperator(tau/2., mu);

  	  // Loop through time steps 
      for(int i = 1; i < m_-1; i++) {
        // Reaction for next full time step 
        sol_next = cap.cwiseQuotient(ones + (cap.cwiseQuotient(sol_cur) - ones) * std::exp(-lambda_*tau));
        sol_cur = sol_next;
        // Diffusion for next full time step 
        sol_next = diffusionEvolutionOperator(tau, sol_cur);
        sol_cur = sol_next;
      }

      // Last half time step: Reaction 
      sol_next = cap.cwiseQuotient(ones + (cap.cwiseQuotient(sol_cur) - ones) * std::exp(-lambda_*tau/2.));
      sol_cur = sol_next;

      sol = sol_cur;
      return sol;
    }


  private:
    // Finite Element Space 
    const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_;
    // Final Time 
    double T_;
    // Number of Timesteps 
    unsigned int m_;
    // Growth Factor 
    double lambda_;
    // coefficient for SDIRK-2 Butcher Tableau
    double xi_;
    // Galerkin matrix corresponding to the negative Laplacian with Robin Boundary Conditions */
    Eigen::SparseMatrix<double> A_;
    // Galerkin matrix for the Mass Matrix 
    Eigen::SparseMatrix<double> M_;
    // Precompute LU decomposition needed for time stepping 
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

};

/* Constructor for StrangSplit */
template <typename DIFF_COEFF>
StrangSplit::StrangSplit(const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space, double T, unsigned m,
                            double lambda, DIFF_COEFF &&c) : fe_space_(fe_space), T_(T), m_(m), lambda_(lambda) {

  const lf::assemble::DofHandler &dofh{fe_space_->LocGlobMap()};

  std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> galerkinpair = assembleGalerkinMatrices(dofh, c);
  A_ = galerkinpair.first;
  M_ = galerkinpair.second;

  // Butcher Tableau Coefficient for SDIRK-2 
  xi_ = 1.0 - 0.5 * sqrt(2.0);

}


} /* namespace FisherKPP */
