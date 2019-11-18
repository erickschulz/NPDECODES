/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

namespace ElementMatrixComputation {

class MyLinearFEElementMatrix {
 public:
  using elem_mat_t = Eigen::Matrix<double, 4, 4>;
  using ElemMat = const elem_mat_t;

  /**
   * @brief Idle constructor
   */
  MyLinearFEElementMatrix() = default;

  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
  /*
   * @brief main routine for the computation of element matrices
   *
   * @param cell reference to the (triangular or quadrilateral) cell for
   *        which the element matrix should be computed.
   * @return a 4x4 matrix, containing the element matrix. The bottom row/column
   *         is not used in the case of a triangle.
   */
  ElemMat Eval(const lf::mesh::Entity &cell);

 private:
  /**
   * @brief LehrFEM's internal laplace element matrix builder
   */
  lf::uscalfe::LinearFELaplaceElementMatrix laplace_elmat_builder_;
};

}  // namespace ElementMatrixComputation
