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

std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<double>>
	CompBeta(std::shared_ptr<const lf::mesh::Mesh> mesh_p, 
			 const Eigen::VectorXd& mu){
  	
	auto beta_p = lf::mesh::utils::make_CodimMeshDataSet(mesh_p, 1, 0.0);

	for(const lf::mesh::Entity* edge: mesh_p->Entities(1)){
		auto endpoints = edge->SubEntities(1);
		unsigned int i = mesh_p->Index(*(endpoints[0]));
		unsigned int j = mesh_p->Index(*(endpoints[1]));

		(*beta_p)(*edge) =  std::exp(mu(i))*Bernoulli(mu(j)-mu(i));
	}

	return beta_p;
}



Eigen::Matrix3d ExpFittedEMP::Eval(const lf::mesh::Entity &cell) {
			
			LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(), "Only 2D triangles are supported.");

			const auto &geom = cell.Geometry();
			double area = lf::geometry::Volume(*geom);
			const auto vertices = lf::geometry::Corners(*geom);

			Eigen::MatrixXd grads(2, 3);
			grads << vertices.col(1) - vertices.col(2),
					 vertices.col(2) - vertices.col(0),
					 vertices.col(0) - vertices.col(1);
			
			Eigen::Matrix3d AK = (1.0 / (4.0 * area)) * grads.transpose() * grads;
			
			Eigen::Vector3d b = beta_loc(cell);

			Eigen::Matrix3d result;
			result << AK(0,1) * b(0) + AK(0,2) * b(2),-AK(0,1) * b(0)                 , -AK(0,2) * b(2),
					 -AK(0,1) * b(0)				 , AK(0,1) * b(0) + AK(1,2) * b(1), -AK(1,2) * b(1),
					 -AK(0,2) * b(2)			     ,-AK(1,2) * b(1)                 ,  AK(0,2) * b(2) + AK(1,2)*b(1);
			
			Eigen::Vector3d mu_exp = (-mu_loc(cell)).array().exp();
			result *= mu_exp.asDiagonal();

			return std::move(result);	
}
  

Eigen::Vector3d ExpFittedEMP::beta_loc(const lf::mesh::Entity& cell){
	Eigen::Vector3d b;
	auto edges = cell.SubEntities(1);
	for(int i = 0; i < 3; ++i){
		b(i) = (*beta_)(*(edges[i]));
	}
	return b;
}

Eigen::Vector3d ExpFittedEMP::mu_loc(const lf::mesh::Entity& cell){
	Eigen::Vector3d m;
	auto mesh_p = fe_space_->Mesh();
	auto vertices = cell.SubEntities(2);
	for(int i = 0; i < 3;++i){
		int index = mesh_p->Index(*(vertices[i]));
		m(i) = mu_(index);
	}
	return m;
}


} /* namespace ExpFittedUpwind */
