/**
 * @file quasiinterpolation.h
 * @brief NPDE exam TEMPLATE HEADER FILE
 * @author Oliver Rietmann
 * @date 15.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef QUASIINTERPOLATION_H_
#define QUASIINTERPOLATION_H_

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include <lf/base/base.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

namespace QuasiInterpolation {

/**
 * @brief Computes the length of the longest edge
 *
 * @param edges range of edges
 * @return length of a longest edge in the range and zero if the range is empty
 */
double maxLength(const nonstd::span<const lf::mesh::Entity *const> &edges);

/**
 * @brief Produces a mapping of vertices $p$ to the corresponding $(K_p,j)$,
 * where $K_p$ is a largest triangle containing $p$ and $j\in\{0,1,2\}$ is the
 * local index of $p$ in $K_p$.
 *
 * @param mesh_p shared pointer to a triangular mesh
 * @return data set storing a mapping as described above
 */
lf::mesh::utils::CodimMeshDataSet<
    std::pair<const lf::mesh::Entity *, unsigned int>>
findKp(std::shared_ptr<const lf::mesh::Mesh> mesh_p);

/**
 * @brief Compute the quasi-interpolation operator $Q_h(v)$.
 *
 * @param fe_space FEM space of tent functions to project onto
 * @param v_mf mesh function representing $v(x)$, such that MESHFUNCTION
 * satisfies the LehrFEM concept MeshFunction of value type R=\texttt{double}
 * @return coefficient vector of size equal to the number of DOFs in fe_space,
 * representing the coefficients of the projected function $Q_h(v)$
 */
/* SAM_LISTING_BEGIN_1 */
template <typename MESHFUNCTION>
Eigen::VectorXd
quasiInterpolate(const lf::uscalfe::FeSpaceLagrangeO1<double> &fe_space,
                 MESHFUNCTION &&v_mf) {
  // Get quadrature points and weights on the reference triangle for a
  // quadrature rule that integrates quadratic polynomials exactly
  lf::quad::QuadRule quadrule =
      lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2);
  const int P = quadrule.NumPoints();
  Eigen::MatrixXd zeta_hat = quadrule.Points();
  Eigen::VectorXd w_hat = quadrule.Weights();
  // ========================================
  // First part of your solution here
  // ========================================
  // Retrieve the map $p \mapsto (K_p, j)$
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space.Mesh();
  lf::mesh::utils::CodimMeshDataSet<
      std::pair<const lf::mesh::Entity *, unsigned int>>
      Kp_mesh_data_set = findKp(mesh_p);

  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};
  const int N = dofh.NumDofs();
  // Vector for returning the basis expansion coefficients of the
  // quasi-interpolant
  Eigen::VectorXd coefficients(N);
  // ========================================
  // Second part of your solution here
  // ========================================
  return coefficients;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Compute the interpolation error for convergence plots given a sequence
 * of $L$ refinements.
 *
 * @param l2_error vector of size $L$ to hold the errors w.r.t. $\lVert\ \cdot\
 * \rVert_{L^2(\Omega)}$.
 * @param h1_error vector of size $L$ to hold the errors w.r.t. $\lVert\ \cdot\
 * \rVert_{H^1(\Omega)}$.
 * @param mesh_hierarchy_p sequence of $L$ refinements of a triangular mesh
 * @param u_mf mesh function representing the function $u(x)$ to project,
 * such that MESHFUNCTION_U satisfies the LehrFEM concept MESHFUNCTION of value
 * type R=\texttt{double}
 * @param grad_u_mf mesh function representing the gradient $\nabla u(x)$ of
 * $u(x)$, such that MESHFUNCTION_GRAD_U satisfies the LehrFEM concept
 * MESHFUNCTION of value type R=\texttt{Eigen::Vector2d}
 */
template <typename MESHFUNCTION_U, typename MESHFUNCTION_GRAD_U>
void interpolationError(
    Eigen::VectorXd &l2_error, Eigen::VectorXd &h1_error,
    std::shared_ptr<lf::refinement::MeshHierarchy> mesh_hierarchy_p,
    MESHFUNCTION_U &&u_mf, MESHFUNCTION_GRAD_U &&grad_u_mf) {
  int L = mesh_hierarchy_p->NumLevels();
  l2_error = Eigen::VectorXd(L);
  h1_error = Eigen::VectorXd(L);

  for (int k = 0; k < L; ++k) {
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = mesh_hierarchy_p->getMesh(k);
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

    Eigen::VectorXd coefficients = quasiInterpolate(*fe_space, u_mf);

    lf::uscalfe::MeshFunctionL2NormDifference loc_comp_l2(fe_space, u_mf, 4);
    l2_error(k) =
        lf::uscalfe::NormOfDifference(fe_space->LocGlobMap(), loc_comp_l2,
                                      coefficients, lf::base::PredicateTrue());

    lf::uscalfe::MeshFunctionL2GradientDifference loc_comp_h1(fe_space,
                                                              grad_u_mf, 4);
    double h1semi_error =
        lf::uscalfe::NormOfDifference(fe_space->LocGlobMap(), loc_comp_h1,
                                      coefficients, lf::base::PredicateTrue());
    auto square = [](double x) { return x * x; };
    h1_error(k) = std::sqrt(square(h1semi_error) + square(l2_error(k)));
  }
}

} // namespace QuasiInterpolation

#endif // #ifndef QUASIINTERPOLATION_H_
