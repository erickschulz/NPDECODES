/**
 * @file minimalgraphsurface.cc
 * @brief NPDE homework 5-3 Minimal Graph Surface code
 * @author R. Hiptmair & W. Tonnonw
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "minimalgraphsurface.h"

#include <lf/base/lf_assert.h>
#include <lf/fe/fe_tools.h>

#include <Eigen/Core>

namespace MinimalGraphSurface {

/* SAM_LISTING_BEGIN_1 */
double computeGraphArea(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    const Eigen::VectorXd& mu_vec) {
  double area;
#if SOLUTION
  // Create a MeshFunction representing the gradient
  lf::fe::MeshFunctionGradFE<double, double> graduh(fes_p, mu_vec);
  // A lambda function realizing a MeshFunction
  // $\sqrt{1+\N{\grad u_h}^2}$
  auto integrand = [&graduh](
                       const lf::mesh::Entity& e,
                       const Eigen::MatrixXd& refc) -> std::vector<double> {
    const std::vector<Eigen::VectorXd> gradvals{graduh(e, refc)};
    std::vector<double> ret(gradvals.size());
    for (int i = 0; i < gradvals.size(); ++i) {
      ret[i] = std::sqrt(1.0 + gradvals[i].squaredNorm());
    }
    return ret;
  };
  area = lf::fe::IntegrateMeshFunction(*fes_p->Mesh(), integrand, 2);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return area;
}
/* SAM_LISTING_END_1 */

// Implementation of the constructor
/* SAM_LISTING_BEGIN_2 */
#if SOLUTION
CoeffTensorA::CoeffTensorA(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    const Eigen::VectorXd& mu)
    : graduh_(fes_p, mu) {}
#else
//====================
// Your code goes here
//====================
#endif
/* SAM_LISTING_END_2 */

// Implementation of evaluation operator for class CoeffTensorA
/* SAM_LISTING_BEGIN_3 */
std::vector<Eigen::Matrix2d> CoeffTensorA::operator()(
    const lf::mesh::Entity& e, const Eigen::MatrixXd& refc) {
  // Number of points for which evaluation is requested
  const int nvals = refc.cols();
  // For returning values
  std::vector<Eigen::Matrix2d> Avals(nvals);
#if SOLUTION
  // Gradients of FE function in those points
  const std::vector<Eigen::VectorXd> gradvals{graduh_(e, refc)};
  LF_ASSERT_MSG_CONSTEXPR(gradvals.size() == nvals,
                          "Wrong number of gradients");
  // Compute tensor A at all input locations
  for (int i = 0; i < nvals; ++i) {
    const Eigen::Vector2d g{gradvals[i]};
    const double norms_g = g.squaredNorm();
    Avals[i] =
        1.0 / (1.0 + norms_g) *
        (Eigen::Matrix2d::Identity() - 2 * g * g.transpose() / (1.0 + norms_g));
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return Avals;
}
/* SAM_LISTING_END_3 */

// Implementation of constructor for CoeffScalarc
/* SAM_LISTING_BEGIN_5 */
#if SOLUTION
CoeffScalarc::CoeffScalarc(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    const Eigen::VectorXd& mu)
    : graduh_(fes_p, mu) {}
#else
//====================
// Your code goes here
//====================
#endif
/* SAM_LISTING_END_5 */

// Implementation of evaluation operator for class CoeffScalarc
/* SAM_LISTING_BEGIN_4 */
std::vector<double> CoeffScalarc::operator()(const lf::mesh::Entity& e,
                                             const Eigen::MatrixXd& refc) {
  // Number of points for which evaluation is requested
  const int nvals = refc.cols();
  // For returning values
  std::vector<double> cvals(nvals);
#if SOLUTION
  // Gradients of FE function in those points
  const std::vector<Eigen::VectorXd> gradvals{graduh_(e, refc)};
  LF_ASSERT_MSG_CONSTEXPR(gradvals.size() == nvals,
                          "Wrong number of gradients");
  // Compute coefficient c for all input points
  for (int i = 0; i < nvals; ++i) {
    cvals[i] = -1.0 / (1.0 + gradvals[i].squaredNorm());
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return cvals;
}
/* SAM_LISTING_END_4 */

}  // namespace MinimalGraphSurface
