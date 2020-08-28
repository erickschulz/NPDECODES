/**
 * @file expfittedupwind.cc
 * @brief NPDE homework ExpFittedUpwind
 * @author Amélie Loher
 * @date 27.08.2020
 * @copyright Developed at ETH Zurich
 */

#include <cmath>
#include <memory>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include "expfittedupwind.h"


namespace ExpFittedUpwind {


double Bernoulli(double tau) {
  
	if(std::abs(tau) < 1e-10) {
		return 1.0;
	} else if(std::abs(tau) < 1e-3) {
		return 1.0 / (1.0 + (0.5 + 1.0/6.0 * tau ) * tau);
	} else {
		return tau / (std::exp(tau) - 1.0);
	}

}

//std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<Eigen::VectorXd>> 
Eigen::VectorXd
	compBeta(const lf::uscalfe::FeSpaceLagrangeO1<double> &fe_space, 
					 const Eigen::VectorXd& mu) {

  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space.Mesh();	
	unsigned int nEdges = mu.size() * (mu.size()-1)/2;
	Eigen::VectorXd beta(nEdges);
  
	unsigned int cnt = 0;
	for(int i = 0; i < mu.size(); i++) {
		for(int j = i+1; j < mu.size(); j++) {
			beta(cnt) = std::exp(mu(i)) * Bernoulli(mu(j) - mu(i));
			++cnt;
		}
	}
  /*
	auto beta_dataset = 
	  lf::mesh::utils::make_CodimMeshDataSet(mesh_p, 1, beta);
  
	return beta_dataset;
	*/
	return beta;
}

class ExpFittedEMP {

	public:
		
		ExpFittedEMP() = delete;
 		explicit ExpFittedEMP(lf::uscalfe::FeSpaceLagrangeO1<double> fe_space,
								 				Eigen::VectorXd mu);
		virtual ~ExpFittedEMP() = default;

		bool isActive(const lf::mesh::Entity &/*cell*/) { return true; }
		
		Eigen::Matrix3d Eval(const lf::mesh::Entity &cell) {
			
			LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(), "Only 2D triangles are supported.");

			const auto &geom = cell.Geometry();
			double area = lf::geometry::Volume(*geom);
			const auto vertices = lf::geometry::Corners(*geom);

			Eigen::MatrixXd grads(2, 3);
			grads << vertices(1) - vertices(2),
							 vertices(2) - vertices(0),
							 vertices(0) - vertices(1);
			Eigen::Matrix3d AK = (1.0 / (4.0 * area)) * grads.transpose() * grads;
			
			auto beta = compBeta(fe_space_, mu_);
			Eigen::Matrix3d result;
			result << AK(0,1) * beta(0) + AK(0,2) * beta(1), -AK(0,1) * beta(0), -AK(0,2) * beta(1),
								-AK(0,1) * beta(0), AK(0,1) * beta(0) + AK(1,2) * beta(2), -AK(1,2) * beta(2),
								-AK(0,2) * beta(1), -AK(1,2) * beta(2), AK(0,2) * beta(1) + AK(1,2) * beta(2);
			
			Eigen::Vector3d mu_exp; mu_exp.array() = mu_.array().exp();
			result *= mu_exp.asDiagonal();

			return std::move(result);	
		}

		private: 
			lf::uscalfe::FeSpaceLagrangeO1<double> fe_space_;
			Eigen::VectorXd mu_;

};
  
template<typename FUNC_F, typename FUNC_G> 
Eigen::VectorXd solveDriftDiffusionDirBVP(const lf::uscalfe::FeSpaceLagrangeO1<double> &fe_space, 
																					const Eigen::VectorXd& mu, FUNC_F &&func_f,FUNC_G &&func_g) {

  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};
	lf::base::size_type N_dofs = dofh.NumDofs();
  
  lf::mesh::utils::MeshFunctionGlobal mf_f{func_f};
	lf::mesh::utils::MeshFunctionGlobal mf_g{func_g};

	// Assembly of the Element Matrix based on our ExpFittedEMP class
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
	ExpFittedEMP elmat_builder(fe_space, mu);
	lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
	
	// Assembly of the RHS
	Eigen::VectorXd phi(N_dofs);
	phi.setZero();
  
	lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)> elvec_builder(fe_space, mf_f);
	lf::assemble::AssembleVectorLocally(0, dofh, elvec_builder, phi);

	// Fixing Solution components according to essential Dirichlet boundary conditions
	auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space.Mesh(), 1)};

	std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<double>> rsf_edge = 
		fe_space.ShapeFunctionLayout(lf::base::RefEl::kSegment());

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


} /* namespace ExpFittedUpwind */
