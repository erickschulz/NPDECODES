/**
 * @file expfittedupwind.h
 * @brief NPDE homework ExpFittedUpwind
 * @author Amélie Loher, Philippe Peter
 * @date 07.01.2021
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <memory>

namespace ExpFittedUpwind {

/**
 * @brief Computes the Bernoulli function B(tau)
 **/
double Bernoulli(double tau);

/**
 * @brief computes the quantities \beta(e) for all the edges e of a mesh
 * @param mesh_p underlying mesh
 * @param mu vector of nodal values of a potential Psi
 * @return  Mesh Data set containing the quantities \beta(e)
 */
std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<double>> CompBeta(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p, const Eigen::VectorXd& mu);

/**
 * @brief computes the element matrices for the exponentially fitted upwind
 * method */
class ExpFittedEMP {
 public:
  /**
   * @brief
   * @param fe_sapce Underlying finite element space
   * @param mu vector of nodal values of a potential Psi
   */
  /* SAM_LISTING_BEGIN_3 */
  explicit ExpFittedEMP(
      std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space,
      Eigen::VectorXd mu)
      : fe_space_(fe_space), beta_(CompBeta(fe_space->Mesh(), mu)), mu_(mu) {}
  /* SAM_LISTING_END_3 */

  virtual ~ExpFittedEMP() = default;

  /**
   * @brief If true, then a triangle is taken into account during assembly
   */
  bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }

  /**
   * @brief actual computation of the element matrix
   * @param cell reference to the triangle for which the matrix is evaluated
   * @return 3x3 dense matrix containg the element matrix
   * */
  Eigen::Matrix3d Eval(const lf::mesh::Entity& cell);

 private:
  /**
   * @brief returns the quanties beta(e) for  the
   * three edges e_0, e_1 and e_2 of a triangle.
   * @param cell reference to the triangle for which the quantities are needed
   * @return vector  [beta(e_0),beta(e_1),beta(e_2)]'
   **/
  Eigen::Vector3d beta_loc(const lf::mesh::Entity& cell);

  /** @brief returns the nodal values of the potential Psi for the
   * three vertices a_1, a_2 and a_3 of a triangle
   * @param cell reference to the triangle for which the quantities are needed
   * @return vector [Psi(a_1), Psi(a_2), Psi(a_3)]'
   **/

  Eigen::Vector3d mu_loc(const lf::mesh::Entity& cell);

  std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>>
      fe_space_;  // underlying finite element space
  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<double>>
      beta_;            // quantities beta(e) for the potential Psi
  Eigen::VectorXd mu_;  // nodal values of the potential Psi
  lf::uscalfe::LinearFELaplaceElementMatrix
      laplace_provider_;  // EMP for A_K, the EM for (u,v) -> \int_\Omega \grad
                          // u * \grad v dx
};

/**
 * @brief computes the approximate solution u_N for the pure Dirichlet BVP
 *        based on the exponentially fitted upwind method.
 *
 * @tparam FUNC_F, FUNC_G  functions that define scalar valued coefficients.
 *
 * @param fe_space underlying finite element space, on which the problem is
 * solved
 * @param mu nodal values of the potential Psi
 * @param func_f source term
 * @param func_g Dirichlet data
 *
 * @return basis expansion coefficient of the approximate solution u_N
 */
template <typename FUNC_F, typename FUNC_G>
Eigen::VectorXd SolveDriftDiffusionDirBVP(
    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    const Eigen::VectorXd& mu, FUNC_F&& func_f, FUNC_G&& func_g) {
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  lf::base::size_type N_dofs = dofh.NumDofs();

  lf::mesh::utils::MeshFunctionGlobal mf_f{func_f};
  lf::mesh::utils::MeshFunctionGlobal mf_g{func_g};

  Eigen::VectorXd sol = Eigen::VectorXd::Ones(N_dofs);

#if SOLUTION
  // Assembly of the System Matrix based on our ExpFittedEMP class
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  ExpFittedEMP elmat_builder(fe_space, mu);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // Assembly of the RHS
  Eigen::VectorXd phi(N_dofs);
  phi.setZero();
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      elvec_builder(fe_space, mf_f);
  lf::assemble::AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // Fixing solution components according to essential Dirichlet boundary
  // conditions
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};

  const lf::fe::ScalarReferenceFiniteElement<double>* rsf_edge =
      fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());

  LF_ASSERT_MSG(rsf_edge != nullptr, "FE specification for edges missing.");

  // Fetch flags and values for dofs on the Dirichlet edges
  auto bd_flags_val{lf::fe::InitEssentialConditionFromFunction(
      *fe_space,
      [&bd_flags](const lf::mesh::Entity& edge) -> bool {
        return bd_flags(edge);
      },
      mf_g)};

  // Eliminate Dirichlet dofs from linear system
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&bd_flags_val](lf::assemble::glb_idx_t idx) {
        return bd_flags_val[idx];
      },
      A, phi);

  // Solve the sparse LSE
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

  solver.compute(A_crs);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  sol = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
#else
  //====================
  // Your code goes here
  //====================
#endif

  return sol;
}

}  // namespace ExpFittedUpwind
