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

#include "preptrimesh.h"

namespace MinimalGraphSurface {

/* SAM_LISTING_BEGIN_1 */
double computeGraphArea(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    const Eigen::VectorXd& mu_vec) {
  double area;
  // Create a MeshFunction representing the gradient
  lf::fe::MeshFunctionGradFE<double, double> graduh(fes_p, mu_vec);
  // A lambda function realizing a MeshFunction $\sqrt{1+\N{\grad u_h}^2}$
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
  const unsigned int quad_ord = 2;
  area = lf::fe::IntegrateMeshFunction(*fes_p->Mesh(), integrand, quad_ord);
  return area;
}
/* SAM_LISTING_END_1 */

// Implementation of the constructor
/* SAM_LISTING_BEGIN_2 */
CoeffTensorA::CoeffTensorA(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    const Eigen::VectorXd& mu)
    : graduh_(fes_p, mu) {}
/* SAM_LISTING_END_2 */

// Implementation of evaluation operator for class CoeffTensorA
/* SAM_LISTING_BEGIN_3 */
std::vector<Eigen::Matrix2d> CoeffTensorA::operator()(
    const lf::mesh::Entity& e, const Eigen::MatrixXd& refc) const {
  // Number of points for which evaluation is requested
  const int nvals = refc.cols();
  // For returning values
  std::vector<Eigen::Matrix2d> Avals(nvals);
  // Gradients of FE function in those points
  const std::vector<Eigen::VectorXd> gradvals{graduh_(e, refc)};
  LF_ASSERT_MSG(gradvals.size() == nvals, "Wrong number of gradients");
  // Compute tensor A at all input locations
  for (int i = 0; i < nvals; ++i) {
    const Eigen::Vector2d g{gradvals[i]};
    const double norms_g = g.squaredNorm();
    Avals[i] =
        1.0 / std::sqrt(1.0 + norms_g) *
        (Eigen::Matrix2d::Identity() - g * g.transpose() / (1.0 + norms_g));
  }
  return Avals;
}
/* SAM_LISTING_END_3 */

// Implementation of constructor for CoeffScalarc
/* SAM_LISTING_BEGIN_5 */
CoeffScalarc::CoeffScalarc(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    const Eigen::VectorXd& mu)
    : graduh_(fes_p, mu) {}
/* SAM_LISTING_END_5 */

// Implementation of evaluation operator for class CoeffScalarc
/* SAM_LISTING_BEGIN_4 */
std::vector<double> CoeffScalarc::operator()(
    const lf::mesh::Entity& e, const Eigen::MatrixXd& refc) const {
  // Number of points for which evaluation is requested
  const int nvals = refc.cols();
  // For returning values
  std::vector<double> cvals(nvals);
  // Gradients of FE function in those points
  const std::vector<Eigen::VectorXd> gradvals{graduh_(e, refc)};
  LF_ASSERT_MSG(gradvals.size() == nvals, "Wrong number of gradients");
  // Compute coefficient c for all input points
  for (int i = 0; i < nvals; ++i) {
    cvals[i] = -1.0 / std::sqrt(1.0 + gradvals[i].squaredNorm());
  }
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
  // Set up an empty sparse matrix to hold the Galerkin matrix
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Initialize ELEMENT_MATRIX_PROVIDER object
  // Tensor coefficient provided by auxiliary MeshFunction object
  CoeffTensorA mf_alpha(fes_p, mu_vec);
  // No zero-order term
  lf::mesh::utils::MeshFunctionConstant<double> mf_zero(0.0);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider elmat_builder(
      fes_p, std::move(mf_alpha), mf_zero);
  // Cell-oriented assembly over the whole computational domain
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
  // Assembly of right-hand-side vector
  Eigen::VectorXd phi(N_dofs);
  {
    // Auxliary sparse matrix
    lf::assemble::COOMatrix<double> T(N_dofs, N_dofs);
    // Scalar coefficient $c(\Bx)$ provided by auxiliary MeshFunction object
    CoeffScalarc mf_c(fes_p, mu_vec);
    lf::uscalfe::ReactionDiffusionElementMatrixProvider Tmat_builder(
        fes_p, std::move(mf_c), mf_zero);
    // Cell-oriented assembly over the whole computational domain
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, Tmat_builder, T);
    phi = T.MatVecMult(1.0, mu_vec);
  }
  // Enforce zero Dirichlet boundary conditions
  // Obtain an array of boolean flags for the edges of the mesh, 'true'
  // indicates that the edge lies on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fes_p->Mesh(), 2)};
  // Elimination of degrees of freedom on the boundary
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&bd_flags,
       &dofh](lf::assemble::glb_idx_t gdof_idx) -> std::pair<bool, double> {
        const lf::mesh::Entity& node{dofh.Entity(gdof_idx)};
        return (bd_flags(node) ? std::make_pair(true, 0.0)
                               : std::make_pair(false, 0.0));
      },
      A, phi);
  // Assembly completed: Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Solve linear system using Eigen's sparse direct elimination
  // Examine return status of solver in case the matrix is singular
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  sol_vec = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
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
      fe_space_p,
      [](Eigen::Vector2d x) -> double { return std::abs(x[0] + x[1] - 1.0); },
      1.0E-4, 1.0E-6, 20,
      [&iterates, &itcnt](const Eigen::VectorXd& mu) -> void {
        iterates.push_back(mu);
        itcnt++;
      }
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
  lf::io::VtkWriter vtk_writer(mesh_p, vtkfile);
  // Create a MeshFunction encoding the FE solution
  lf::fe::MeshFunctionFE<double, double> mf_uh(fe_space_p, mu);
  // Write nodal data to file
  vtk_writer.WritePointData("mPoint", mf_uh);
  std::cout << "\n The finite element solution was written to:" << std::endl;
  std::cout << ">> " << vtkfile << std::endl;
  // Supplement: output for visualization with MATLAB
  PrepTriMesh::prepTriMesh(fe_space_p, mu, "minsurf");
  /* MATLAB code for visualization:
     [x,y,TRI] = minsurf();
     nodvals_minsurf;
     trimesh(TRI(:,1:3),x,y,nodevals);
     xlabel("{\bf x}");
     ylabel("{\bf y}");
     zlabel("{\bf u}");
     print -depsc2 'minsurf.eps';
   */
}
/* SAM_LISTING_END_7 */

}  // namespace MinimalGraphSurface
