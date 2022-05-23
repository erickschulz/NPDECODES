/**
 * @file transpsemilagr.h
 * @brief NPDE homework TranspSemiLagr code
 * @author Tobias Rohner
 * @date November 2020
 * @copyright Developed at SAM, ETH Zurich
 */
#include "transpsemilagr.h"

namespace TranspSemiLagr {

void enforce_zero_boundary_conditions(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    lf::assemble::COOMatrix<double>& A, Eigen::VectorXd& b) {
  lf::mesh::utils::MeshFunctionGlobal mf_zero{
      [](const Eigen::Vector2d& /*x*/) { return 0.0; }};
  const lf::fe::ScalarReferenceFiniteElement<double>* rsf_edge_p =
      fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());

  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
  auto flag_values{
      lf::fe::InitEssentialConditionFromFunction(*fe_space, bd_flags, mf_zero)};

  lf::assemble::FixFlaggedSolutionCompAlt<double>(
      [&flag_values](lf::assemble::glb_idx_t dof_idx) {
        return flag_values[dof_idx];
      },
      A, b);
}

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd solverot(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    Eigen::VectorXd u0_vector, int N, double T) {
  //====================
  // Your code goes here
  //====================
  return (T + N) * Eigen::VectorXd::Ones(u0_vector.size());
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solvetrp(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    Eigen::VectorXd u0_vector, int N, double T) {
  //====================
  // Your code goes here
  //====================
  return (T + N) * Eigen::VectorXd::Ones(u0_vector.size());
}
/* SAM_LISTING_END_2 */

}  // namespace TranspSemiLagr
