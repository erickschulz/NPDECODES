/**
 * @file stableevaluationatapoint.h
 * @brief NPDE homework StableEvaluationAtAPoint
 * @author Amélie Loher
 * @date 22/04/2020
 * @copyright Developed at ETH Zurich
 */

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <cmath>
#include <complex>

namespace StableEvaluationAtAPoint {

/* Returns the mesh size for the given mesh. */
double getMeshSize(const std::shared_ptr<const lf::mesh::Mesh> &mesh_p) {
  double mesh_size = 0.0;

  // Find maximal edge length
  double edge_length;
  for (const lf::mesh::Entity *edge : mesh_p->Entities(1)) {
    // Compute the length of the edge
    auto endpoints = lf::geometry::Corners(*(edge->Geometry()));
    edge_length = (endpoints.col(0) - endpoints.col(1)).norm();
    if (mesh_size < edge_length) {
      mesh_size = edge_length;
    }
  }

  return mesh_size;
}

/* Returns G(x,y). */
double G(Eigen::Vector2d x, Eigen::Vector2d y) {
  double res;

  LF_ASSERT_MSG(x != y, "G not defined for these coordinates!");
  res = (-1.0 / (2.0 * M_PI)) * std::log((x - y).norm());

  return res;
}

/* Returns the gradient of G(x,y). */
Eigen::Vector2d gradG(Eigen::Vector2d x, Eigen::Vector2d y) {
  Eigen::Vector2d res;

  LF_ASSERT_MSG(x != y, "G not defined for these coordinates!");

  res = (x - y) / (2.0 * M_PI * (x - y).squaredNorm());

  return res;
}

/* Evaluates the Integral P_SL using the local midpoint rule
 * on the partitioning of the boundary of Omega induced by the mesh.
 * The supplied meshes are unitary squares.
 */
template <typename FUNCTOR>
double PSL(std::shared_ptr<const lf::mesh::Mesh> mesh, const FUNCTOR &v,
           const Eigen::Vector2d x) {
  double PSLval = 0.0;

  // This predicate returns true for edges on the boundary
  auto bd_flags_edge{lf::mesh::utils::flagEntitiesOnBoundary(mesh, 1)};

  // Loop over boundary edges
  for (const lf::mesh::Entity *e : mesh->Entities(1)) {
    if (bd_flags_edge(*e)) {
      const lf::geometry::Geometry *geo_ptr = e->Geometry();
      LF_ASSERT_MSG(geo_ptr != nullptr, "Missing geometry!");
      // Fetch coordinates of corner points
      Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
      // Determine midpoints of edges
      Eigen::Vector2d midpoint;
      midpoint(0) = 0.5 * (corners(0, 0) + corners(0, 1));
      midpoint(1) = 0.5 * (corners(1, 0) + corners(1, 1));

      // Compute and add the elemental contribution
      PSLval += v(midpoint) * G(x, midpoint) * lf::geometry::Volume(*geo_ptr);
    }
  }

  return PSLval;
}

/* Evaluates the Integral P_DL using the local midpoint rule
 * on the partitioning of the boundary of Omega induced by the mesh.
 * The supplied meshes are unitary squares.
 */
template <typename FUNCTOR>
double PDL(std::shared_ptr<const lf::mesh::Mesh> mesh, const FUNCTOR &v,
           const Eigen::Vector2d x) {
  double PDLval = 0.0;

  // This predicate returns true for edges on the boundary
  auto bd_flags_edge{lf::mesh::utils::flagEntitiesOnBoundary(mesh, 1)};

  // Loop over boundary edges
  for (const lf::mesh::Entity *e : mesh->Entities(1)) {
    if (bd_flags_edge(*e)) {
      const lf::geometry::Geometry *geo_ptr = e->Geometry();
      LF_ASSERT_MSG(geo_ptr != nullptr, "Missing geometry!");
      // Fetch coordinates of corner points
      Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
      // Determine midpoints of edges
      Eigen::Vector2d midpoint;
      midpoint(0) = 0.5 * (corners(0, 0) + corners(0, 1));
      midpoint(1) = 0.5 * (corners(1, 0) + corners(1, 1));
      // Determine the normal vector n on the unit square
      Eigen::Vector2d n;
      if ((midpoint(0) > midpoint(1)) && (midpoint(0) < 1.0 - midpoint(1))) {
        n << 0.0, -1.0;

      } else if ((midpoint(0) > midpoint(1)) &&
                 (midpoint(0) > 1.0 - midpoint(1))) {
        n << 1.0, 0.0;

      } else if ((midpoint(0) < midpoint(1)) &&
                 (midpoint(0) > 1.0 - midpoint(1))) {
        n << 0.0, 1.0;

      } else {
        n << -1.0, 0.0;
      }

      // Compute and add the elemental contribution
      PDLval += v(midpoint) * (gradG(x, midpoint)).dot(n) *
                lf::geometry::Volume(*geo_ptr);
    }
  }

  return PDLval;
}

/* This function computes u(x) = P_SL(grad u * n) - P_DL(u).
 * For u(x) = log( (x + (1, 0)^T).norm() ) and x = (0.3, 0.4)^T,
 * it computes the difference between the analytical and numerical
 * evaluation of u. The mesh is supposed to be the unit square.
 */
double pointEval(std::shared_ptr<const lf::mesh::Mesh> mesh) {
  double error = 0.0;

  const auto u = [](Eigen::Vector2d x) -> double {
    Eigen::Vector2d one(1.0, 0.0);
    return std::log((x + one).norm());
  };

  const auto gradu = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    Eigen::Vector2d one(1.0, 0.0);
    return (x + one) / (x + one).squaredNorm();
  };

  // Define a Functor for the dot product of grad u(x) * n(x)
  const auto dotgradu_n = [gradu](const Eigen::Vector2d x) -> double {
    // Determine the normal vector n on the unit square
    Eigen::Vector2d n;
    if ((x(0) > x(1)) && (x(0) < 1.0 - x(1))) {
      n << 0.0, -1.0;

    } else if ((x(0) > x(1)) && (x(0) > 1.0 - x(1))) {
      n << 1.0, 0.0;

    } else if ((x(0) < x(1)) && (x(0) > 1.0 - x(1))) {
      n << 0.0, 1.0;

    } else {
      n << -1.0, 0.0;
    }

    return gradu(x).dot(n);
  };

  // Compute right hand side
  const Eigen::Vector2d x(0.3, 0.4);
  const double rhs = PSL(mesh, dotgradu_n, x) - PDL(mesh, u, x);
  // Compute the error
  error = std::abs(u(x) - rhs);

  return error;
}

/* Computes Psi_x(y). */
double Psi(const Eigen::Vector2d y) {
  double Psi_xy;

  Eigen::Vector2d half(0.5, 0.5);
  double constant = M_PI / (0.5 * std::sqrt(2) - 1.0);
  double dist = (y - half).norm();

  if (dist <= 0.25 * std::sqrt(2)) {
    Psi_xy = 0.0;

  } else if (dist >= 0.5) {
    Psi_xy = 1.0;

  } else {
    Psi_xy = std::pow(std::cos(constant * (dist - 0.5)), 2);
  }

  return Psi_xy;
}

/* Computes grad(Psi_x(y)). */
Eigen::Vector2d gradPsi(const Eigen::Vector2d y) {
  Eigen::Vector2d gradPsi_xy;

  Eigen::Vector2d half(0.5, 0.5);
  double constant = M_PI / (0.5 * std::sqrt(2) - 1.0);
  double dist = (y - half).norm();

  if (dist <= 0.25 * std::sqrt(2)) {
    gradPsi_xy(0) = 0.0;
    gradPsi_xy(1) = 0.0;

  } else if (dist >= 0.5) {
    gradPsi_xy(0) = 0.0;
    gradPsi_xy(1) = 0.0;

  } else {
    gradPsi_xy = -2.0 * std::cos(constant * (dist - 0.5)) *
                 std::sin(constant * (dist - 0.5)) * (constant / dist) *
                 (y - half);
  }

  return gradPsi_xy;
}

/* Computes Laplacian of Psi_x(y). */
double laplPsi(const Eigen::Vector2d y) {
  double laplPsi_xy;

  Eigen::Vector2d half(0.5, 0.5);
  double constant = M_PI / (0.5 * std::sqrt(2) - 1.0);
  double dist = (y - half).norm();

  if (dist <= 0.25 * std::sqrt(2)) {
    laplPsi_xy = 0.0;

  } else if (dist >= 0.5) {
    laplPsi_xy = 0.0;

  } else {
    laplPsi_xy =
        (2 * std::pow(constant, 2) / (y - half).squaredNorm()) *
            (y - half).dot(y - half) *
            (std::pow(std::sin(constant * (dist - 0.5)), 2) -
             std::pow(std::cos(constant * (dist - 0.5)), 2)) -
        (2 * constant / dist) * std::cos(constant * (dist - 0.5)) *
            std::sin(constant * (dist - 0.5)) *
            (1.0 - ((y - half).dot(y - half) / (y - half).squaredNorm()));
  }

  return laplPsi_xy;
}

/* Computes Jstar
 * fe_space: finite element space defined on a triangular mesh of the square
 * domain u: Function handle for u x: Coordinate vector for x
 */
template <typename FUNCTOR>
double Jstar(std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space,
             FUNCTOR &&u, const Eigen::Vector2d x) {
  double val = 0.0;

  std::shared_ptr<const lf::mesh::Mesh> mesh = fe_space->Mesh();

  // Use midpoint quadrature rule
  const lf::quad::QuadRule qr = lf::quad::make_TriaQR_MidpointRule();
  // Quadrature points
  const Eigen::MatrixXd zeta_ref{qr.Points()};
  // Quadrature weights
  const Eigen::VectorXd w_ref{qr.Weights()};
  // Number of quadrature points
  const lf::base::size_type P = qr.NumPoints();

  // Loop over all cells
  for (const lf::mesh::Entity *entity : mesh->Entities(0)) {
    LF_ASSERT_MSG(entity->RefEl() == lf::base::RefEl::kTria(),
                  "Not on triangular mesh!");

    const lf::geometry::Geometry &geo{*entity->Geometry()};
    const Eigen::MatrixXd zeta{geo.Global(zeta_ref)};
    const Eigen::VectorXd gram_dets{geo.IntegrationElement(zeta_ref)};

    for (int l = 0; l < P; l++) {
      val -= w_ref[l] * u(zeta.col(l)) *
             (2.0 * (gradG(x, zeta.col(l))).dot(gradPsi(zeta.col(l))) +
              G(x, zeta.col(l)) * laplPsi(zeta.col(l))) *
             gram_dets[l];
    }
  }

  /* VARIANT:
  auto lambda = lf::mesh::utils::MeshFunctionGlobal( [&] (Eigen::Vector2d y) {

        return (-u(y) * (2.0 * gradG(x, y).dot(gradPsi(y)) + G(x, y) *
  laplPsi(y) ));
  }
  );

  double val_test = lf::uscalfe::IntegrateMeshFunction(*mesh, lambda, 9);
  */
  return val;
}

/* Evaluates u(x) according to (3.11.14).
 * u: Function Handle for u
 * x: Coordinate vector for x
 */
template <typename FUNCTOR>
double stab_pointEval(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space,
    FUNCTOR &&u, const Eigen::Vector2d x) {
  double res = 0.0;

  Eigen::Vector2d half(0.5, 0.5);
  if ((x - half).norm() <= 0.25) {
    res = Jstar(fe_space, u, x);

  } else {
    std::cerr << "The point does not fulfill the assumptions" << std::endl;
  }

  return res;
}

} /* namespace StableEvaluationAtAPoint */
