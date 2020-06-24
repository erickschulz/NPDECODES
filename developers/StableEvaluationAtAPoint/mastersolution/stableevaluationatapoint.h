/**
 * @file stableevaluationatapoint.h
 * @brief NPDE homework StableEvaluationAtAPoint
 * @author Amélie Loher & Erick Schulz
 * @date 22/04/2020
 * @copyright Developed at ETH Zurich
 */

#include <cmath>
#include <complex>

#include <Eigen/Core>

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>

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
  // Straightforward implementation
  res = (-1.0 / (2.0 * M_PI)) * std::log((x - y).norm());
  return res;
}

/* Returns the gradient of G(x,y). */
Eigen::Vector2d gradG(Eigen::Vector2d x, Eigen::Vector2d y) {

  Eigen::Vector2d res;
  LF_ASSERT_MSG(x != y, "G not defined for these coordinates!");
  // Straightforward implementation
  res = (x - y) / (2.0 * M_PI * (x - y).squaredNorm());
  return res;
}

/* Evaluates the Integral P_SL using the local midpoint rule
 * on the partitioning of the boundary of Omega induced by the mesh.
 * The supplied meshes are unitary squares.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR>
double PSL(std::shared_ptr<const lf::mesh::Mesh> mesh, FUNCTOR &&v,
           const Eigen::Vector2d x) {
  double PSLval = 0.0;
#if SOLUTION
  // This predicate returns true for edges on the boundary
  auto bd_flags_edge{lf::mesh::utils::flagEntitiesOnBoundary(mesh, 1)};
  // Loop over boundary edges
  for (const lf::mesh::Entity *e : mesh->Entities(1)) {
    if (bd_flags_edge(*e)) {
      const lf::geometry::Geometry *geo_ptr = e->Geometry();
      LF_ASSERT_MSG(geo_ptr != nullptr, "Missing geometry!");
      // Fetch coordinates of corner points
      const Eigen::Matrix2d corners = lf::geometry::Corners(*geo_ptr);
      // Determine midpoints of edges
      const Eigen::Vector2d midpoint{0.5 * (corners.col(0) + corners.col(1))};
      // Compute and add the elemental contribution
      PSLval += v(midpoint) * G(x, midpoint) * lf::geometry::Volume(*geo_ptr);
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return PSLval;
}

/* SAM_LISTING_END_1 */

/* Evaluates the Integral P_DL using the local midpoint rule
 * on the partitioning of the boundary of Omega induced by the mesh.
 * The supplied meshes are unitary squares.
 */
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
double PDL(std::shared_ptr<const lf::mesh::Mesh> mesh, FUNCTOR &&v,
           const Eigen::Vector2d x) {
  double PDLval = 0.0;
#if SOLUTION
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
      // Determine the normal vector n on the unit square. There are four
      // possibilities depending on the edge of the unit square
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
#else
  //====================
  // Your code goes here
  //====================
#endif
  return PDLval;
}

/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
/* This function computes u(x) = P_SL(grad u * n) - P_DL(u).
 * For u(x) = log( (x + (1, 0)^T).norm() ) and x = (0.3, 0.4)^T,
 * it computes the difference between the analytical and numerical
 * evaluation of u. The mesh is supposed to cover the unit square.
 */
double pointEval(std::shared_ptr<const lf::mesh::Mesh> mesh) {
  double error = 0.0;
#if SOLUTION
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
#else
  //====================
  // Your code goes here
  //====================
#endif
  return error;
}

/* SAM_LISTING_END_3 */

/* Computes Psi_x(y). */
double Psi(const Eigen::Vector2d y) {
  double Psi_xy;
  const Eigen::Vector2d half(0.5, 0.5);
  const double constant = M_PI / (0.5 * std::sqrt(2) - 1.0);
  const double dist = (y - half).norm();

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
/* SAM_LISTING_BEGIN_4 */
double Jstar(std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
             Eigen::VectorXd uFE, const Eigen::Vector2d x) {
  double val = 0.0;

#if SOLUTION
  std::shared_ptr<const lf::mesh::Mesh> mesh = fe_space->Mesh();

  // Use midpoint quadrature rule
  const lf::quad::QuadRule qr = lf::quad::make_TriaQR_MidpointRule();
  // Quadrature points
  const Eigen::MatrixXd zeta_ref{qr.Points()};
  // Quadrature weights
  const Eigen::VectorXd w_ref{qr.Weights()};
  // Number of quadrature points
  const lf::base::size_type P = qr.NumPoints();

  // Create mesh function to be evaluated at the quadrature points
  auto uFE_mf = lf::uscalfe::MeshFunctionFE(fe_space, uFE);

  // Loop over all cells
  for (const lf::mesh::Entity *entity : mesh->Entities(0)) {
    LF_ASSERT_MSG(entity->RefEl() == lf::base::RefEl::kTria(),
                  "Not on triangular mesh!");

    const lf::geometry::Geometry &geo{*entity->Geometry()};
    const Eigen::MatrixXd zeta{geo.Global(zeta_ref)};
    const Eigen::VectorXd gram_dets{geo.IntegrationElement(zeta_ref)};

    auto u_vals = uFE_mf(*entity, zeta_ref);

    for (int l = 0; l < P; l++) {
      val -= w_ref[l] * u_vals[l] *
             (2.0 * (gradG(x, zeta.col(l))).dot(gradPsi(zeta.col(l))) +
              G(x, zeta.col(l)) * laplPsi(zeta.col(l))) *
             gram_dets[l];
    }
  }

#else
  //====================
  // Your code goes here
  //====================
#endif
  return val;
}

/* SAM_LISTING_END_4 */

/* Evaluates u(x) according to (3.11.14).
 * u: Function Handle for u
 * x: Coordinate vector for x
 */
double
stab_pointEval(std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
               Eigen::VectorXd uFE, const Eigen::Vector2d x) {

  double res = 0.0;

#if SOLUTION
  Eigen::Vector2d half(0.5, 0.5);
  if ((x - half).norm() <= 0.25) {
    res = Jstar(fe_space, uFE, x);

  } else {
    std::cerr << "The point does not fulfill the assumptions" << std::endl;
  }

#else
//====================
// Your code goes here
//====================
#endif
  return res;
}

template <typename FUNCTOR>
Eigen::VectorXd solveBVP(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p,
    FUNCTOR &&u) {
  Eigen::VectorXd discrete_solution;

  // TOOLS AND DATA
  // Pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  // Obtain specification for shape functions on edges
  std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<double>>
      rsf_edge_p = fe_space_p->ShapeFunctionLayout(lf::base::RefEl::kSegment());

  // Dirichlet data
  auto mf_g = lf::mesh::utils::MeshFunctionGlobal(
      [u](Eigen::Vector2d x) -> double { return u(x); });
  // Right-hand side source function f
  auto mf_f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  // I : ASSEMBLY
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Right hand side vector, must be initialized with 0!
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();

  // I.i : Computing volume matrix for negative Laplace operator
  // Initialize object taking care of local mass (volume) computations.
  lf::uscalfe::LinearFELaplaceElementMatrix elmat_builder{};
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrix A.
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // I.ii : Computing right-hand side vector
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      elvec_builder(fe_space_p, mf_f);
  // Invoke assembly on cells (codim == 0)
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // I.iii : Imposing essential boundary conditions
  // Obtain an array of boolean flags for the edges of the mesh, 'true'
  // indicates that the edge lies on the boundary (codim = 1)
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Inspired by the example in the documentation of
  // InitEssentialConditionFromFunction()
  // https://craffael.github.io/lehrfempp/namespacelf_1_1uscalfe.html#a5afbd94919f0382cf3fb200c452797ac
  // Determine the fixed dofs on the boundary and their values
  // Alternative: See lecturedemoDirichlet() in
  // https://github.com/craffael/lehrfempp/blob/master/examples/lecturedemos/lecturedemoassemble.cc
  auto edges_flag_values_Dirichlet{
      lf::uscalfe::InitEssentialConditionFromFunction(dofh, *rsf_edge_p,
                                                      bd_flags, mf_g)};
  // Eliminate Dirichlet dofs from the linear system
  lf::assemble::FixFlaggedSolutionCompAlt<double>(
      [&edges_flag_values_Dirichlet](lf::assemble::glb_idx_t gdof_idx) {
        return edges_flag_values_Dirichlet[gdof_idx];
      },
      A, phi);

  // Assembly completed! Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_sparse = A.makeSparse();

  // II : SOLVING  THE LINEAR SYSTEM
  // II.i : Setting up Eigen's sparse direct elimination
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_sparse);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  // II.ii : Solving
  discrete_solution = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

  return discrete_solution;
}; // solveBVP

} /* namespace StableEvaluationAtAPoint */
