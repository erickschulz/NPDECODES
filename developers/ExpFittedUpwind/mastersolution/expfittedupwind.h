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

//std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<Eigen::VectorXd>>
Eigen::VectorXd
	compBeta(std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> &fe_space,
           const Eigen::VectorXd& mu);

template<typename FUNC_F, typename FUNC_G>
Eigen::VectorXd solveDriftDiffusionDirBVP(std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> &fe_space,
                                          const Eigen::VectorXd& mu, FUNC_F &&func_f,FUNC_G &&func_g);

}
