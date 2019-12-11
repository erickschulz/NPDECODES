/**
 * @file
 * @brief NPDE homework SimpleLinearFiniteElements
 * @author Amélie Loher
 * @date 11/12/2019
 * @copyright Developed at ETH Zurich
 */

#include <functional>
#include <string>

#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "tria_mesh_2D.h"

namespace SimpleLinearFiniteElements {

Eigen::Matrix<double, 2, 3> gradbarycoordinates(const Eigen::Matrix<double, 2, 3>& Vertices);

Eigen::Matrix3d ElementMatrix_Mass_LFE(const Eigen::Matrix<double, 2, 3>& Vertices);

Eigen::Matrix3d ElementMatrix_Lapl_LFE(const Eigen::Matrix<double, 2, 3>& Vertices);

Eigen::Matrix3d ElementMatrix_LaplMass_LFE(const Eigen::Matrix<double, 2, 3>& Vertices);

Eigen::Vector3d localLoadLFE(const Eigen::Matrix<double, 2, 3>& Vertices,
                             const std::function<double(const Eigen::Vector2d&)>& FHandle);

double L2Error(const SimpleLinearFiniteElements::TriaMesh2D& mesh, const Eigen::VectorXd& uFEM,
               const std::function<double(double, double)> exact);

double H1Serror(const SimpleLinearFiniteElements::TriaMesh2D& mesh, const Eigen::VectorXd& uFEM,
                const std::function<Eigen::Vector2d(double, double)> exact);

Eigen::SparseMatrix<double> GalerkinAssembly(
    const SimpleLinearFiniteElements::TriaMesh2D& Mesh, const std::function<Eigen::Matrix3d(const Eigen::Matrix<double, 2, 3>&)>& getElementMatrix);

Eigen::VectorXd assemLoad_LFE(const SimpleLinearFiniteElements::TriaMesh2D& Mesh,
                              const std::function<Eigen::Vector3d(const Eigen::Matrix<double, 2, 3>&, std::function<double(const Eigen::Vector2d&)>)>& getElementVector,
                              const std::function<double(const Eigen::Vector2d&)>& FHandle);


} // namespace SimpleLinearFiniteElements 
