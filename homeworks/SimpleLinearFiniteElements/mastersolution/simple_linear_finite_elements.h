#include <functional>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace SimpleLinearFiniteElements {

using TriGeo_t = Eigen::Matrix<double, 2, 3>;
using FHandle_t = std::function<double(const Eigen::Vector2d&)>;
using LocalMatrixHandle_t = std::function<Eigen::Matrix3d(const TriGeo_t&)>;
using LocalVectorHandle_t =
    std::function<Eigen::Vector3d(const TriGeo_t&, FHandle_t)>;
using Triplet_t = std::vector<Eigen::Triplet<double>>;

Eigen::Matrix3d ElementMatrix_Mass_LFE(const TriGeo_t& vertices);

/**
 * @brief simple mesh data structure used for this problem
 */
struct TriaMesh2D {
  // Constructor: reads mesh data from file
  TriaMesh2D(const std::string&);
  virtual ~TriaMesh2D(void) {}

  // Creates EPS rendering of mesh geometry using MathGL
  void plotMesh(const std::string& epsfile, int drawvertices = 0) const;
  void plotSurf(const std::string& epsfile,
                const Eigen::VectorXd& values) const;

  // Data members describing geometry and topolgy
  Eigen::Matrix<double, Eigen::Dynamic, 2> Coordinates;
  Eigen::Matrix<int, Eigen::Dynamic, 3> Elements;
};

double L2Error(const TriaMesh2D& mesh, const Eigen::VectorXd& uFEM,
               const std::function<double(double, double)> exact);

double H1Serror(const TriaMesh2D& mesh, const Eigen::VectorXd& uFEM,
                const std::function<Eigen::Vector2d(double, double)> exact);

Eigen::Matrix<double, 2, 3> gradbarycoordinates(const TriGeo_t& Vertices);

Eigen::Vector3d localLoadLFE(const TriGeo_t& Vertices,
                             const FHandle_t& FHandle);

Eigen::Matrix3d ElementMatrix_Lapl_LFE(const TriGeo_t& Vertices);

Eigen::Matrix3d ElementMatrix_LaplMass_LFE(const TriGeo_t& Vertices);

Eigen::SparseMatrix<double> GalerkinAssembly(
    const TriaMesh2D& Mesh, const LocalMatrixHandle_t& getElementMatrix);

Eigen::VectorXd assemLoad_LFE(const TriaMesh2D& Mesh,
                              const LocalVectorHandle_t& getElementVector,
                              const FHandle_t& FHandle);

}  // namespace SimpleLinearFiniteElements
