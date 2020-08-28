/**
 * @file expfittedupwind_main.cc
 * @brief NPDE homework ExpFittedUpwind
 * @author Amélie Loher
 * @date 27.08.2020
 * @copyright Developed at ETH Zurich
 */

#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <lf/mesh/hybrid2d/hybrid2d.h>

#include "expfittedupwind.h"

int main() {
  // Obtain Mesh on unit square
	auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
	const lf::io::GmshReader reader(std::move(mesh_factory),
                                  CURRENT_SOURCE_DIR "/../meshes/square.msh");
	auto mesh_p = reader.mesh();
	// FE Space, Dofhandler
	auto fe_space = std::make_shared<FeSpaceLagrangeO2<double>>(mesh_p);
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  
	// Output file
	std::ofstream L2output;
  L2output.open("L2error.txt");
  L2output << "No. of dofs, L2 error" << std::endl;
  
	// Data
	Eigen::Vector2d q = Eigen::Vector2d::Ones(2);
	auto Psi = [&q] (Eigen::Vector2d x) -> double {
		return q.dot(x);
	};

	auto f = [] (Eigen::Vector2d x) -> double { return 0.0; };

	auto g = [&Psi] (Eigen::Vector2d x) -> double {
		return Psi(x).exp();
	};
  
  Eigen::VectorXd mu(N_dofs);
	mu.setZero();

	for(const lf::mesh::Entity *el: mesh.Entities(0)) {
		const auto& geom = el.Geometry();
		const auto ids = dofh.GlobalDofIndices(*el);
		const lf::assemble::size_type N_locDofs(dofh.NumLocalDofs(*el));
		for(const auto idx: ids) {

			Eigen::MatrixXd vertices(2,3);
			auto vertices = lf::geometry::Corners(*geom);
			for(int loc_idx = 0; loc_idx < N_locDofs; loc_idx++) {
				Eigen::Vector2d pt = vertices.col(loc_idx);
			}
			mu(idx) = pt.x() + pt.y();
		}
	}

  // Compute Solution
	Eigen::VectorXd sol_vec = ExpFittedUpwind::solveDriftDiffusionBVP(fe_space, 
																																mu, f, g);
  
	// auto sol = MeshFunctionFE(fe_space, sol_vec);

  // Compare to exact solution interpolated on our mesh
  auto ref_sol = [&Psi] (Eigen::Vector2d x) -> double {
		return Psi(x).exp();
	};
	lf::mesh::utils::MeshFunctionGlobal mf_ref_sol{ref_sol};
	Eigen::VectorXd ref_sol_vec(N_dofs);
	ref_sol_vec = lf::uscalfe::NodalProjection(fe_space, mf_ref_sol);
  
  double L2_error = (sol_vec - ref_sol_vec).norm() / N_dofs;
  L2output << N_dofs << ", " << L2_error << std::endl;

  return 0;
}
