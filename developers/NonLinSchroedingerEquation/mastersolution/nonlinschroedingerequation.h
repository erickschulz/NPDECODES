/**
 * @file nonlinschroedingerequation.h
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 22.04.2020
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <lf/mesh/mesh.h>

namespace NonLinSchroedingerEquation {

class MassElementMatrixProvider {
public:
  /** @brief Default implement: all cells are active */
  bool isActive(const lf::mesh::Entity &cell) { return true; }
  /*
   * @brief routine for the computation of element matrices
   * @param cell reference to the triangular cell for
   *        which the element matrix should be computed.
   * @return element matrix
   */
  Eigen::Matrix3d Eval(const lf::mesh::Entity &cell);
};

class Energy {
public:
  Energy() { }
  Energy(const Eigen::SparseMatrix<double> &A, const Eigen::SparseMatrix<double> &D) : _A(A), _D(D) { }

  double operator()(const Eigen::VectorXcd &mu) const;
private:
  Eigen::SparseMatrix<double> _A;
  Eigen::SparseMatrix<double> _D;
};

/*class StiffnessElementMatrixProvider {
public:
  bool isActive(const lf::mesh::Entity &cell) { return true; }
  Eigen::Matrix3d Eval(const lf::mesh::Entity &cell);
};

class TrapezoidalElementVectorProvider {
public:
  bool isActive(const lf::mesh::Entity &cell) { return true; }
  Eigen::Vector3d Eval(const lf::mesh::Entity &cell);
};*/

}  // namespace NonLinSchroedingerEquation
