#ifndef CD_TOOLS_H
#define CD_TOOLS_H

/**
 * @file cd_tools.h
 * @brief Utility Functions for the Convection-Diffusion Problem
 * @author Philippe Peter
 * @date July 2021
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace ConvectionDiffusion {

/**
 *@brief Computes the diameter of a TRIANGULAR mesh entity
 */
double Diameter(const lf::mesh::Entity& entity);

/**
 *@brief Computes the mesh width of a TRIANGULAR mesh
 */
double MeshWidth(std::shared_ptr<const lf::mesh::Mesh> mesh_p);

/**
 * @brief Evaluates a MeshFunction at a point specified by its global
 * coordinates
 */
template <typename MF>
double EvaluateMeshFunction(std::shared_ptr<const lf::mesh::Mesh> mesh_p, MF mf,
                            Eigen::Vector2d global, double tol = 10E-10) {
  for (const lf::mesh::Entity* entity_p : mesh_p->Entities(0)) {
    LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity_p->RefEl(),
                  "Function only defined for triangular cells");

    // compute geometric information about the cell
    const lf::geometry::Geometry* geo_p = entity_p->Geometry();
    Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

    // transform global coordinates to local coordinates on the cell
    Eigen::Matrix2d A;
    A << corners.col(1) - corners.col(0), corners.col(2) - corners.col(0);
    Eigen::Vector2d b;
    b << global - corners.col(0);
    Eigen::Vector2d loc = A.fullPivLu().solve(b);

    // evaluate meshfunction, if local coordinates lie in the reference triangle
    if (loc(0) >= 0 - tol && loc(1) >= 0 - tol && loc(0) + loc(1) <= 1 + tol) {
      return mf(*entity_p, loc)[0];
    }
  }
  return 0.0;
}
/**
 * @brief Evaluates a MeshFunction at points specified by their global
 * coordinates
 * @param mesh_p a TRIANGULAR mesh
 */
template <typename MF>
std::vector<double> EvaluateMeshFunction(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p, MF mf,
    const std::vector<Eigen::Vector2d>& global, double tol = 10E-10) {
  unsigned N = global.size();
  std::vector<double> res(N);
  std::vector<bool> computed(N, false);
  unsigned counter = 0;

  for (const lf::mesh::Entity* entity_p : mesh_p->Entities(0)) {
    // verify that entity_p is a triangle
    LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity_p->RefEl(),
                  "Function only defined for triangular cells");

    // compute geometric information about the cell
    const lf::geometry::Geometry* geo_p = entity_p->Geometry();
    Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

    // Prepare global to local coordinate transformation
    Eigen::Matrix2d A;
    A << corners.col(1) - corners.col(0), corners.col(2) - corners.col(0);
    Eigen::FullPivLU<Eigen::Matrix2d> lu = A.fullPivLu();

    for (unsigned i = 0; i < N; ++i) {
      if (computed[i] == false) {
        // compute local coordinates
        Eigen::Vector2d b = global[i] - corners.col(0);
        Eigen::Vector2d loc = lu.solve(b);

        // evaluate MeshFunction, if local coordinates lie in the reference
        // triangle
        if (loc(0) >= 0 - tol && loc(1) >= 0 - tol &&
            loc(0) + loc(1) <= 1 + tol) {
          res[i] = mf(*entity_p, loc)[0];
          computed[i] = true;
          counter++;

          // Check, if MeshFunction was evaluated at all points
          if (counter == N) {
            return res;
          }
        }
      }
    }
  }
  return res;
}

/**
 * @brief Samples a MeshFunction at points along a curve
 * @param file_name Outut file
 * @param mesh_p underlying mesh
 * @param gamma Curve parametrized over the unit interval [0,1]
 * @param mf MeshFunction
 * @param N Number of sample points
 */
template <typename CURVE, typename MF>
void SampleMeshFunction(std::string file_name,
                        std::shared_ptr<const lf::mesh::Mesh> mesh_p,
                        CURVE gamma, MF mf, int N) {
  // Sample uniformly along gamma
  Eigen::VectorXd sample_times = Eigen::VectorXd::LinSpaced(N, 0.0, 1.0);
  std::vector<Eigen::Vector2d> sample_points(N);
  for (int i = 0; i < N; ++i) {
    double t = sample_times(i);
    sample_points[i] = gamma(t);
  }

  std::vector<double> res = EvaluateMeshFunction(mesh_p, mf, sample_points);

  // output
  std::ofstream file;
  file.open(file_name);
  for (int i = 0; i < N; ++i) {
    double t = sample_times(i);
    Eigen::Vector2d x = sample_points[i];
    double eval = res[i];
    file << t << ", " << x(0) << ", " << x(1) << ", " << eval << std::endl;
  }
  file.close();
}

}  // namespace ConvectionDiffusion

#endif  // CD_TOOLS_H