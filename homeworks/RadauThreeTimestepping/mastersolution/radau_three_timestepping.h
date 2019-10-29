#ifndef RADAU_H
#define RADAU_H

/** @file
 * @brief NPDE RadauThreeTimestepping
 * @author Erick Schulz
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */
#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <boost/assert.hpp>
#include <boost/filesystem.hpp>

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <unsupported/Eigen/KroneckerProduct>

namespace RadauThreeTimestepping {

/** @brief time depedent heat source */
Eigen::VectorXd rhsVectorheatSource(const lf::assemble::DofHandler &,
                                    double time);

/** @brief solve heat equation with rhsVectorHeat source as source */
Eigen::VectorXd solveHeatEvolution(const lf::assemble::DofHandler &dofh,
                                   unsigned int m, double final_time);

/** @Brief This function enforces Dirichlet zero boundary conditions on the
 * Galerkin stiffness and mass matrices. It transforms every columns and vectors
 * associated to a global index belonging to a degree of freedom lying on the
 * boundary to zero entries but the diagonal one which is set to 1.0
 * @param selectvals is the predicate identifying the boundary indices of the
 * rows and columns that are to be dropped */
template <typename SCALAR, typename SELECTOR>
void dropMatrixRowsColumns(SELECTOR &&selectvals,
                           lf::assemble::COOMatrix<SCALAR> &A) {
  const lf::assemble::size_type N(A.cols());
  LF_ASSERT_MSG(A.rows() == N, "Matrix must be square!");
  A.setZero(
      [&selectvals](lf::assemble::gdof_idx_t i, lf::assemble::gdof_idx_t j) {
        return (selectvals(i) || selectvals(j));
      });
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N; ++dofnum) {
    const auto selval{selectvals(dofnum)};
    if (selval) {
      A.AddToEntry(dofnum, dofnum, 1.0);
    }
  }
}

/** @brief This class implements a Lehrfem++ matrix provider defining a
LinFEMassMatrixProvider::Eval function returning the local MASS matrix for
linear first-order lagrange FE bases over triangular mesh (only!). Integration
over the triangular cells is performed using the trapezoidal rule.*/
class LinFEMassMatrixProvider {
 public:
  /** @brief default constructor */
  explicit LinFEMassMatrixProvider() = default;
  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
  /** @brief Main method for computing the element vector
   * @param cell refers to current cell for which the element vector is desired
   * The implementation uses an analytic formula defined over triangular cells*/
  Eigen::Matrix<double, 3, 3> Eval(const lf::mesh::Entity &tria);
};  // class LinFEMassMatrixProvider

/** @brief This class implements a Lehrfem++ matrix provider defining a
TrapRuleLinFEElemVecProvider<FUNCTOR>::Eval function returning the local
contribution to the element vectors for linear first-order lagrange FE bases
over triangular mesh (only!). Integration over the triangular cells is performed
using the trapezoidal rule.*/
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>  // lambda predicate
class TrapRuleLinFEElemVecProvider {
 public:
  /** @brief Constructor storing the right hand side function */
  explicit TrapRuleLinFEElemVecProvider(FUNCTOR f) : f_(f) {}
  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
  /** @brief Main method for computing the element vector
   * @param cell current cell for which the element vector is desired
   * The implementation uses simple vertex based quadrature and an approximation
   * of the volume of a cell just using the integration element at the
   * barycenter.*/
  Eigen::Vector3d Eval(const lf::mesh::Entity &tria);

 private:
  /* f_ provides the evaluation of the source function at coordinates*/
  FUNCTOR f_;
};
/* SAM_LISTING_END_2 */

// TrapRuleLinFEElemVecProvider
/* Implementing member function Eval of class TrapRuleLinFEElemVecProvider*/
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR>
Eigen::Vector3d TrapRuleLinFEElemVecProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity &tria) {
  Eigen::Vector3d ElemVec;
  /* SOLUTION_BEGIN */
  // Throw error in case no triangular cell
  LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << tria.RefEl());
  // Obtain vertex coordinates of the triangle in a 2x3 matrix
  const auto corners{lf::geometry::Corners(*(tria.Geometry()))};
  // Compute the scaling factor for the local load vector
  const double area_third = lf::geometry::Volume(*(tria.Geometry())) / 3.0;
  LF_ASSERT_MSG((corners.cols() == 3) && (corners.rows() == 2),
                "Invalid vertex coordinate " << corners.rows() << "x"
                                             << corners.cols() << " matrix");
  ElemVec = Eigen::Vector3d(area_third * f_(corners.col(0)),
                            area_third * f_(corners.col(1)),
                            area_third * f_(corners.col(2)));
  /* SOLUTION_END */
  return ElemVec;
}  // TrapRuleLinFEElemVecProvider<FUNCTOR>::Eval
/* SAM_LISTING_END_3 */

/** @brief class providing timestepping for heat equation */
class Radau3MOLTimestepper {
 public:
  // Disabled constructors
  Radau3MOLTimestepper() = delete;
  Radau3MOLTimestepper(const Radau3MOLTimestepper &) = delete;
  Radau3MOLTimestepper(Radau3MOLTimestepper &&) = delete;
  Radau3MOLTimestepper &operator=(const Radau3MOLTimestepper &) = delete;
  Radau3MOLTimestepper &operator=(const Radau3MOLTimestepper &&) = delete;
  // Main constructor; precomputations are done here
  Radau3MOLTimestepper(const lf::assemble::DofHandler &dofh);

  // Destructor
  virtual ~Radau3MOLTimestepper() = default;

  /* Class member functions */
  // Discrete evolution operator for Radau IIA 3rd order
  Eigen::VectorXd discreteEvolutionOperator(double time, double tau,
                                            const Eigen::VectorXd &mu) const;

 private:
  /* SOLUTION_BEGIN */
  // const shared_ptr(lf::assemble::DofHandler) dofh_;
  double tau_;
  const lf::assemble::DofHandler &dofh_;  // dangerous
  // Matrices in triplet format holding Galerkin matrices
  Eigen::SparseMatrix<double> A_;     // Element matrix
  Eigen::SparseMatrix<double> A_Kp_;  // Element Kronecker product matrix
  Eigen::SparseMatrix<double> M_Kp_;  // Mass Kronecker product matrix
  // Butcher tableau of the Runge-Kutta RADAU-2 method
  Eigen::Matrix<double, 2, 2> U_;
  Eigen::Vector2d b_;
  Eigen::Vector2d c_;
  // For fixed step-size in time, the linear system of equations implicitely
  // defining the Runge-Kutta increments is independent of time. We can thus
  // precompute the LU decomposition for more efficiency.
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_;
  /* SOLUTION_END */
};  // class Radau3MOLTimestepper

}  // namespace RadauThreeTimestepping

#endif
