/**
 * @file quasiinterpolation_main.cc
 * @brief NPDE exam TEMPLATE MAIN FILE
 * @author Oliver Rietmann
 * @date 15.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <cmath>

#include "iohelper.h"
#include "quasiinterpolation.h"

constexpr double PI = 3.14159265358979323846;
// A bump function
double u1(Eigen::Vector2d x) {
  double abs_x = x.norm();
  return abs_x < 0.5 ? std::cos(PI * abs_x) : 0.0;
};
// The gradient of the bump function
Eigen::Vector2d grad_u1(Eigen::Vector2d x) {
  double abs_x = x.norm();
  if (abs_x == 0.0 || abs_x >= 0.5)
    return Eigen::Vector2d::Zero();
  else
    return -PI * std::sin(PI * abs_x) / abs_x * x;
};

// A simplex auxiliary function
inline double square(double y) { return y * y; }

// A smoother bump function
double u2(Eigen::Vector2d x) { return square(u1(x)); };
// ... and its gradient
Eigen::Vector2d grad_u2(Eigen::Vector2d x) { return 2.0 * u1(x) * grad_u1(x); };

int main(int /* argc */, char** /*argv*/) {
  // Build mesh hierarchy
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory),
                            CURRENT_SOURCE_DIR "/../meshes/mesh.msh");
  std::shared_ptr<lf::mesh::Mesh> base_mesh_p{reader.mesh()};
  std::shared_ptr<lf::refinement::MeshHierarchy> mesh_hierarchy_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(base_mesh_p, 6);

  // Compute meshwidth for each refinement
  int L = mesh_hierarchy_p->NumLevels();
  Eigen::VectorXd meshwidth(L);
  for (int k = 0; k < L; ++k) {
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = mesh_hierarchy_p->getMesh(k);
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    meshwidth(k) = QuasiInterpolation::maxLength(mesh_p->Entities(1));
  }
  {
    // Errors for u1
    lf::mesh::utils::MeshFunctionGlobal u1_mf(u1);
    lf::mesh::utils::MeshFunctionGlobal grad_u1_mf(grad_u1);
    // Retrieve error norms
    auto [l2_error, h1_error] = QuasiInterpolation::interpolationError(
        mesh_hierarchy_p, u1_mf, grad_u1_mf);
    QuasiInterpolation::printError(meshwidth, l2_error, h1_error,
                                   "errors for u1");
    QuasiInterpolation::writeCSV(meshwidth, l2_error, h1_error,
                                 CURRENT_BINARY_DIR "/convergence_u1.csv");
    // Call a Python script to generate plots
    std::system("python3 " CURRENT_SOURCE_DIR
                "/plot_convergence.py " CURRENT_BINARY_DIR
                "/convergence_u1.csv " CURRENT_BINARY_DIR
                "/convergence_u1.eps 'Interpolation error for $u(x)=u_1(x)$'");
  }
  {
    // Errors for u2
    lf::mesh::utils::MeshFunctionGlobal u2_mf(u2);
    lf::mesh::utils::MeshFunctionGlobal grad_u2_mf(grad_u2);
    auto [l2_error, h1_error] = QuasiInterpolation::interpolationError(
        mesh_hierarchy_p, u2_mf, grad_u2_mf);
    QuasiInterpolation::printError(meshwidth, l2_error, h1_error,
                                   "errors for u2");
    QuasiInterpolation::writeCSV(meshwidth, l2_error, h1_error,
                                 CURRENT_BINARY_DIR "/convergence_u2.csv");
    std::system("python3 " CURRENT_SOURCE_DIR
                "/plot_convergence.py " CURRENT_BINARY_DIR
                "/convergence_u2.csv " CURRENT_BINARY_DIR
                "/convergence_u2.eps 'Interpolation error for $u(x)=u_2(x)$'");
  }
  return 0;
}
