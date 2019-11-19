/**
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss
 * @date   16.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "cr_reference_finite_element.h"

namespace NonConformingCrouzeixRaviartFiniteElements {

/* SAM_LISTING_BEGIN_1 */
lf::base::RefEl CRReferenceFiniteElement::RefEl() const {
  lf::base::RefElType ref_el_type;

  //====================
  // Your code goes here
  // TODO: task 2-14.q)
  //====================

  return lf::base::RefEl(ref_el_type);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
unsigned int CRReferenceFiniteElement::Degree() const {
  unsigned int degree;

  //====================
  // Your code goes here
  // TODO: task 2-14.q)
  //====================


  return degree;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
size_type CRReferenceFiniteElement::NumRefShapeFunctions() const {
  size_type num_ref_shape_functions;

  //====================
  // Your code goes here
  // TODO: task 2-14.q)
  //====================

  return num_ref_shape_functions;
}

size_type CRReferenceFiniteElement::NumRefShapeFunctions(dim_t codim) const {
  switch (codim) {
  //====================
  // Your code goes here
  // TODO: task 2-14.q)
  //====================
    default:
      LF_VERIFY_MSG(false, "Codimension out of range for triangle")
      return 0;
  }
}

size_type CRReferenceFiniteElement::NumRefShapeFunctions(
    dim_t codim, sub_idx_t subidx) const {

  switch (codim) {
  //====================
  // Your code goes here
  // TODO: task 2-14.q)
  //====================
    default:
      LF_VERIFY_MSG(false, "Codimension out of range for triangle")
      return 0;
  }
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>
CRReferenceFiniteElement::EvalReferenceShapeFunctions(
    const Eigen::MatrixXd& refcoords) const {
  const auto num_points = static_cast<size_type>(refcoords.cols());
  Eigen::MatrixXd eval_ref_shape_functions(3, num_points);

  //====================
  // Your code goes here
  // TODO: task 2-14.q)
  //====================

  return eval_ref_shape_functions;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>
CRReferenceFiniteElement::GradientsReferenceShapeFunctions(
    const Eigen::MatrixXd& refcoords) const {
  const auto num_points = static_cast<size_type>(refcoords.cols());
  Eigen::MatrixXd grad_ref_shape_functions(3, 2 * num_points);

  //====================
  // Your code goes here
  // TODO: task 2-14.r)
  //====================

  return grad_ref_shape_functions;
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
Eigen::MatrixXd CRReferenceFiniteElement::EvaluationNodes() const {
  Eigen::MatrixXd eval_nodes(2, 3);

  //====================
  // Your code goes here
  // TODO: task 2-14.s)
  //====================

  return eval_nodes;
}
/* SAM_LISTING_END_6 */

size_type CRReferenceFiniteElement::NumEvaluationNodes() const { return 3; }

/* SAM_LISTING_BEGIN_7 */
Eigen::Matrix<scalar_type, 1, Eigen::Dynamic>
CRReferenceFiniteElement::NodalValuesToDofs(
    const Eigen::Matrix<scalar_type, 1, Eigen::Dynamic>& nodvals) const {
  LF_VERIFY_MSG(nodvals.cols() == NumEvaluationNodes(),
                "nodvals = " << nodvals << " <-> " << NumEvaluationNodes());

  Eigen::MatrixXd coeffs;

  //====================
  // Your code goes here
  // TODO: task 2-14.s)
  //====================

  return coeffs;
}
/* SAM_LISTING_END_7 */

}  // namespace NonConformingCrouzeixRaviartFiniteElements
