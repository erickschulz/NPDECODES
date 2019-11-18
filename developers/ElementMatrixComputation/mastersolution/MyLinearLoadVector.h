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

class MyLinearLoadVector {
 public:
  using FHandle_t = std::function<double(const Eigen::Vector2d &)>;
  using elem_vec_t = Eigen::Matrix<double, 4, 1>;
  using ElemVec = const elem_vec_t;

  /** @brief Constructor storing the right hand side function */
  explicit MyLinearLoadVector(FHandle_t f) : f_(f) {}

  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
  /**
   * @brief Main method for computing the element vector
   *
   * @param cell current cell for which the element vector is desired
   *
   * The implementation uses simple edge midpoint based quadrature and an
   * approximation of the volume of a cell just using the integration element at
   * the barycenter.
   */
  ElemVec Eval(const lf::mesh::Entity &cell);

 private:
  /** `f_(x)` where `x` is a 2D vector that provides the evaluation of the
   * source function */
  FHandle_t f_;
};

MyLinearLoadVector::ElemVec computeLoadVector(
    Eigen::MatrixXd vertices, const MyLinearLoadVector::FHandle_t f);

}  // namespace ElementMatrixComputation
