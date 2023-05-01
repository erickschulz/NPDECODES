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
#include <lf/fe/mesh_function_fe.h>

#include <Eigen/Core>


namespace MinimalGraphSurface {

/* SAM_LISTING_BEGIN_1 */
double computeGraphArea(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    const Eigen::VectorXd& mu_vec) {
  double area;
  //====================
  // Your code goes here
  //====================
  return area;
}
/* SAM_LISTING_END_1 */

// Implementation of the constructor
/* SAM_LISTING_BEGIN_2 */
CoeffTensorA::CoeffTensorA(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    const Eigen::VectorXd& mu)
    : graduh_(fes_p, mu) {
  //====================
  // Your code goes here
  // or even elsewhere or nowhere
  //====================
}
/* SAM_LISTING_END_2 */

// Implementation of evaluation operator for class CoeffTensorA
/* SAM_LISTING_BEGIN_3 */
std::vector<Eigen::Matrix2d> CoeffTensorA::operator()(
    const lf::mesh::Entity& e, const Eigen::MatrixXd& refc) const {
  // Number of points for which evaluation is requested
  const int nvals = refc.cols();
  // For returning values
  std::vector<Eigen::Matrix2d> Avals(nvals);
  //====================
  // Your code goes here
  //====================
  return Avals;
}
/* SAM_LISTING_END_3 */

// Implementation of constructor for CoeffScalarc
/* SAM_LISTING_BEGIN_5 */
CoeffScalarc::CoeffScalarc(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    const Eigen::VectorXd& mu)
    : graduh_(fes_p, mu) {
  //====================
  // Your code goes here
  // or somewhere else or nowhere
  //====================
}
/* SAM_LISTING_END_5 */

// Implementation of evaluation operator for class CoeffScalarc
/* SAM_LISTING_BEGIN_4 */
std::vector<double> CoeffScalarc::operator()(
    const lf::mesh::Entity& e, const Eigen::MatrixXd& refc) const {
  // Number of points for which evaluation is requested
  const int nvals = refc.cols();
  // For returning values
  std::vector<double> cvals(nvals);
  //====================
  // Your code goes here
  //====================
  return cvals;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_6 */
Eigen::VectorXd computeNewtonCorrection(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    const Eigen::VectorXd& mu_vec) {
  // Obtain reference to the underlying finite element mesh
  const lf::mesh::Mesh& mesh{*fes_p->Mesh()};
  // The local-to-global index mapping
  const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
  // Get the number of degrees of freedom = dimension of FE space
  const lf::base::size_type N_dofs(dofh.NumDofs());
  LF_ASSERT_MSG(mu_vec.size() == N_dofs, "Vector length mismatch!");
  // Solution vector = return value
  Eigen::VectorXd sol_vec(N_dofs);
  //====================
  // Your code goes here
  //====================
  return sol_vec;
}
/* SAM_LISTING_END_6 */

/* SAM_LISTING_BEGIN_7 */
void graphMinSurfVis(std::string meshfile, std::string vtkfile) {
  // Read mesh for unit square from file
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), meshfile.c_str());
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
  // Finite element space
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Solve non-linear BVP by means of Newton's method
  std::vector<Eigen::VectorXd> iterates{};
  int itcnt = 0;
  Eigen::VectorXd mu = graphMinimalSurface(
      //====================
      // Replace this line with meaningful function argument
      //====================
      fe_space_p, [](Eigen::Vector2d x) { return 0.0; }, 0.0, 0.0, 0
  );
  // Tabulate progress of iteration
  unsigned int N = mu.size();
  for (int i = 0; i < itcnt; ++i) {
    double mu_norm = iterates[i].norm() / std::sqrt(N);
    std::cout << "k = " << i << ": |mu(" << i << ")| = " << mu_norm;
    if (i > 0) {
      std::cout << ", |correction| = "
                << (iterates[i] - iterates[i - 1]).norm() / std::sqrt(N);
    }
    std::cout << std::endl;
  }
  // Output solution for visualization
  //====================
  // Your code goes here
  //====================
}
/* SAM_LISTING_END_7 */

}  // namespace MinimalGraphSurface
