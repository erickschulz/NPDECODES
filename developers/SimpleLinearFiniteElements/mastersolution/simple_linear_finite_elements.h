#include <functional>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "tria_mesh_2D.h"

namespace SimpleLinearFiniteElements {

using TriGeo_t = Eigen::Matrix<double, 2, 3>;
using FHandle_t = std::function<double(const Eigen::Vector2d&)>;
using LocalMatrixHandle_t = std::function<Eigen::Matrix3d(const TriGeo_t&)>;
using LocalVectorHandle_t =
    std::function<Eigen::Vector3d(const TriGeo_t&, FHandle_t)>;
using Triplet_t = std::vector<Eigen::Triplet<double>>;

Eigen::Matrix3d ElementMatrix_Mass_LFE(const TriGeo_t& vertices);

double L2Error(const SimpleLinearFiniteElements::TriaMesh2D& mesh, const Eigen::VectorXd& uFEM,
               const std::function<double(double, double)> exact);

double H1Serror(const SimpleLinearFiniteElements::TriaMesh2D& mesh, const Eigen::VectorXd& uFEM,
                const std::function<Eigen::Vector2d(double, double)> exact);

//Eigen::Matrix<double, 2, 3> gradbarycoordinates(const TriGeo_t& Vertices);

//Eigen::Vector3d localLoadLFE(const TriGeo_t& Vertices,
  //                           const FHandle_t& FHandle);

//Eigen::Matrix3d ElementMatrix_Lapl_LFE(const TriGeo_t& Vertices);

//Eigen::Matrix3d ElementMatrix_LaplMass_LFE(const TriGeo_t& Vertices);

Eigen::SparseMatrix<double> GalerkinAssembly(
    const SimpleLinearFiniteElements::TriaMesh2D& Mesh, const LocalMatrixHandle_t& getElementMatrix);

Eigen::VectorXd assemLoad_LFE(const SimpleLinearFiniteElements::TriaMesh2D& Mesh,
                              const LocalVectorHandle_t& getElementVector,
                              const FHandle_t& FHandle);

}  // namespace SimpleLinearFiniteElements
