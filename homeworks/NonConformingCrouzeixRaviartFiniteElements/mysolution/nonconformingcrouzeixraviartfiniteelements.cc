/** @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 *  @author Anian Ruoss, edited Amélie Loher
 *  @date   16.03.2019, 03.03.20
 *  @copyright Developed at ETH Zurich */

#include "nonconformingcrouzeixraviartfiniteelements.h"

namespace NonConformingCrouzeixRaviartFiniteElements {

// Crouzeix-Raviart finite element space defined on triangular meshes only
lf::base::RefEl CRReferenceFiniteElement::RefEl() const {
  lf::base::RefElType ref_el_type;
  // TODO: task 2-14.q)
  //====================
  // Your code goes here
  //====================
  return lf::base::RefEl(ref_el_type);
}
// Crouzeix-Raviart are piecewise linear polynomials
unsigned int CRReferenceFiniteElement::Degree() const { 
	unsigned int degree;
	// TODO: task 2-14.q)
  	//====================
  	// Your code goes here
  	//====================
	return degree;
}

lf::assemble::size_type CRReferenceFiniteElement::NumRefShapeFunctions() const {
  lf::assemble::size_type num_ref_shape_functions;
  // TODO: task 2-14.q)
    //====================
    // Your code goes here
    //====================
  return num_ref_shape_functions;
}

lf::assemble::size_type CRReferenceFiniteElement::NumRefShapeFunctions(lf::assemble::dim_t codim) const {
  switch (codim) {
	// TODO: task 2-14.q)
    //====================
    // Your code goes here
    //====================
    default:
      LF_VERIFY_MSG(false, "Codimension out of range for triangle")
      return 0;
  }
}

lf::assemble::size_type CRReferenceFiniteElement::NumRefShapeFunctions(
    lf::assemble::dim_t codim, lf::base::sub_idx_t subidx) const {
  switch (codim) {
	// TODO: task 2-14.q)
    //====================
    // Your code goes here
    //====================
    default:
      LF_VERIFY_MSG(false, "Codimension out of range for triangle")
      return 0;
  }
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
CRReferenceFiniteElement::EvalReferenceShapeFunctions(
    const Eigen::MatrixXd& refcoords) const {
  // Data
  const auto num_points = static_cast<lf::assemble::size_type>(refcoords.cols());
  // Initialize a matrix that will store the values of the reference basis
  // functions evaluated at the coordinates passed as arguments
  Eigen::MatrixXd eval_ref_shape_functions(3, num_points);
  // TODO: task 2-14.q)
    //====================
    // Your code goes here
    //====================
  return eval_ref_shape_functions;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
CRReferenceFiniteElement::GradientsReferenceShapeFunctions(
    const Eigen::MatrixXd& refcoords) const {
  // Data
  const auto num_points = static_cast<lf::assemble::size_type>(refcoords.cols());
  // Initialize a matrix that will store the gradients of the reference basis
  // functions evaluated at the coordinates passed as arguments
  Eigen::MatrixXd grad_ref_shape_functions(3, 2 * num_points);
  // TODO: task 2-14.r)
  //====================
  // Your code goes here
  //====================
  return grad_ref_shape_functions;
}

Eigen::MatrixXd CRReferenceFiniteElement::EvaluationNodes() const {
  Eigen::MatrixXd eval_nodes(2, 3);
  // TODO: task 2-14.s)
  //====================
  // Your code goes here
  //====================
  return eval_nodes;
}

lf::assemble::size_type CRReferenceFiniteElement::NumEvaluationNodes() const { 
	return 3; 
}

Eigen::Matrix<double, 1, Eigen::Dynamic>
CRReferenceFiniteElement::NodalValuesToDofs(
    const Eigen::Matrix<double, 1, Eigen::Dynamic>& nodvals) const {
  LF_VERIFY_MSG(nodvals.cols() == NumEvaluationNodes(),
                "nodvals = " << nodvals << " <-> " << NumEvaluationNodes());

  Eigen::MatrixXd coeffs;
  // TODO: task 2-14.s)
  //====================
  // Your code goes here
  //====================
  return nodvals;
}

}  // namespace NonConformingCrouzeixRaviartFiniteElements
