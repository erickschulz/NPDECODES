/**
 * @file MinimalGrapdSurface.h
 * @brief NPDE homework 5-3 Minimal Graph Surface code
 * @author Minimal Graph Surface
 * @date April 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <string>

namespace MinimalGraphSurface {
/** @brief HW 5-3f) Computes the area of a graph
 *
 * @param fes_p pointer to linear Lagrangian finite element space
 * @param mu_vec basis expansion coefficient vector for FE function
 *
 * The area formula is
 * \f[
 * A = \int\limits_{\Omega} \sqrt{1+\left\|\mathfb{grad}\,u_{h}\right\|^2}
 * \f]
 * It can be computed approximated by mesh-based local integration using some
 * local quadrature rule.
 */
double computeGraphArea(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    const Eigen::VectorXd& mu_vec);

/** @brief Tensor coefficient for linearized minimal surface veriational
 * equation:
 *
 * See sub-problem 5-3e)
 */
class CoeffTensorA {
 public:
  CoeffTensorA(
      std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
      const Eigen::VectorXd& mu);
  CoeffTensorA() = delete;
  CoeffTensorA(const CoeffTensorA&) = default;
  CoeffTensorA(CoeffTensorA&&) = default;
  CoeffTensorA& operator=(const CoeffTensorA&) = delete;
  std::vector<Eigen::Matrix2d> operator()(const lf::mesh::Entity& e,
                                          const Eigen::MatrixXd& refc) const;

 private:
  lf::fe::MeshFunctionGradFE<double, double> graduh_;
};

/** @brief Tensor coefficient for linearized minimal surface veriational
 * equation:
 *
 * See sub-problem 5-3e)
 */
class CoeffScalarc {
 public:
  CoeffScalarc(
      std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
      const Eigen::VectorXd& mu);
  CoeffScalarc() = delete;
  CoeffScalarc(const CoeffScalarc&) = default;
  CoeffScalarc(CoeffScalarc&&) = default;
  CoeffScalarc& operator=(const CoeffScalarc&) = delete;
  std::vector<double> operator()(const lf::mesh::Entity& e,
                                 const Eigen::MatrixXd& refc) const;

 private:
  lf::fe::MeshFunctionGradFE<double, double> graduh_;
};

/** @brief Solves discrete variational equation for Newton correction
 *
 * Sub-problem 5-3h)
 */
Eigen::VectorXd computeNewtonCorrection(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    const Eigen::VectorXd& mu_vec);

/** @brief Computation of function whose graph has minimal area over
 * \f$\Omega\f$
 *
 * Sub-problem 5-3i)
 */
/* SAM_LISTING_BEGIN_9 */
template <typename D_FUNCTOR,
          typename RECORDER = std ::function<void(const Eigen::VectorXd&)>>
Eigen::VectorXd graphMinimalSurface(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    D_FUNCTOR&& boundary_data, double rtol, double atol, unsigned int itmax,
    RECORDER&& rec = [](const Eigen::VectorXd&) -> void {}) {
  // Obtain reference to the underlying finite element mesh
  const lf::mesh::Mesh& mesh{*fes_p->Mesh()};
  // The local-to-global index mapping
  const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
  // Get the number of degrees of freedom = dimension of FE space
  const lf::base::size_type N_dofs(dofh.NumDofs());

  // I. Solve Dirichlet problem for $-\Delta u=0$ in order to obtain a good
  // initial guess
  Eigen::VectorXd uh(N_dofs);
  //====================
  // Your code goes here
  //====================
  // II. Newton iteration with correction based termination
  //====================
  // Your code goes here
  //====================
  return uh;
}
/* SAM_LISTING_END_9 */

/** @brief VTK Output of graph with minimal area
 *
 * Sub-problem 5-3j)
 */
void graphMinSurfVis(std::string meshfile, std::string vtkfile);

}  // namespace MinimalGraphSurface
