/**
 * @file expfittedupwind.h
 * @brief NPDE homework ExpFittedUpwind
 * @author Amélie Loher
 * @date 27.08.2020
 * @copyright Developed at ETH Zurich
 */

#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <lf/uscalfe/uscalfe.h>

namespace ExpFittedUpwind {

double Bernoulli(double tau);

std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<double>>
	CompBeta(std::shared_ptr<const lf::mesh::Mesh> mesh_p,
              const Eigen::VectorXd& mu);


class ExpFittedEMP {
 public:
		
 		explicit ExpFittedEMP(std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space,
                              std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<double>> beta,
							  Eigen::VectorXd mu): fe_space_(fe_space),beta_(beta), mu_(mu)  { }
		
		virtual ~ExpFittedEMP() = default;

		bool isActive(const lf::mesh::Entity &/*cell*/) { return true; }
		
        Eigen::Matrix3d Eval(const lf::mesh::Entity& cell); 

private: 
    Eigen::Vector3d beta_loc(const lf::mesh::Entity& cell);

    Eigen::Vector3d mu_loc(const lf::mesh::Entity& cell);

    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_;
    std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<double>> beta_;
    Eigen::VectorXd mu_;

};


template<typename FUNC_F, typename FUNC_G> 
Eigen::VectorXd SolveDriftDiffusionDirBVP(std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space, 
										const Eigen::VectorXd& mu, FUNC_F &&func_f,FUNC_G &&func_g) {

    const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
	lf::base::size_type N_dofs = dofh.NumDofs();
  
   lf::mesh::utils::MeshFunctionGlobal mf_f{func_f};
	lf::mesh::utils::MeshFunctionGlobal mf_g{func_g};

	// Assembly of the Element Matrix based on our ExpFittedEMP class
    lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
    auto beta = CompBeta(fe_space->Mesh(),mu);
	ExpFittedEMP elmat_builder(fe_space, beta, mu);
	lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
	
	// Assembly of the RHS
	Eigen::VectorXd phi(N_dofs);
	phi.setZero();
  
	lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)> elvec_builder(fe_space, mf_f);
	lf::assemble::AssembleVectorLocally(0, dofh, elvec_builder, phi);

	// Fixing Solution components according to essential Dirichlet boundary conditions
	auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};

	std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<double>> rsf_edge = 
		fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());

	LF_ASSERT_MSG(rsf_edge != nullptr, "FE specification for edges missing.");
  
	// Fetch flags and values for dofs on Dirichlet edges
	auto bd_flags_val{lf::uscalfe::InitEssentialConditionFromFunction(dofh, *rsf_edge,
		[&bd_flags](const lf::mesh::Entity& edge) -> bool {
			return bd_flags(edge);
		}, mf_g)};

	// Eliminate Dirichlet dofs from linear system
  lf::assemble::FixFlaggedSolutionComponents<double>(
		[&bd_flags_val](lf::assemble::glb_idx_t idx) {
			return bd_flags_val[idx];
		},
		A, phi);

   // Solve the Sparse LSE
   Eigen::SparseMatrix<double> A_crs = A.makeSparse();
   Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  
	solver.compute(A_crs);
	LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
	Eigen::VectorXd sol = solver.solve(phi);
	LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

    return sol;
}
} //namespace ExpFittedUpwind