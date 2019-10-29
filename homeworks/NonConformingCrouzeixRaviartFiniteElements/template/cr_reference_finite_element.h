/**
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss
 * @date   16.03.2019
 * @copyright Developed at ETH Zurich
 */

#ifndef NUMPDE_CRREFERENCEFINITEELEMENT_H
#define NUMPDE_CRREFERENCEFINITEELEMENT_H

#include <lf/base/base.h>
#include <lf/uscalfe/uscalfe.h>
#include "cr_types.h"

namespace NonConformingCrouzeixRaviartFiniteElements {

class CRReferenceFiniteElement final
    : public lf::uscalfe::ScalarReferenceFiniteElement<scalar_type> {
 public:
  /**
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::RefEl()
   */
  lf::base::RefEl RefEl() const override;

  /**
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::Degree()
   */
  unsigned int Degree() const override;

  /**
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  size_type NumRefShapeFunctions() const override;

  /**
   * @copydoc
   * lf::uscalfe::ScalarReferenceFiniteElement::NumRefShapeFunctions(codim)
   */
  size_type NumRefShapeFunctions(dim_t codim) const override;

  /**
   * @copydoc
   * lf::uscalfe::ScalarReferenceFiniteElement::NumRefShapeFunctions(codim,
   * subidx)
   */
  size_type NumRefShapeFunctions(dim_t codim, sub_idx_t subidx) const override;

  /**
   * @copydoc
   * lf::uscalfe::ScalarReferenceFiniteElement::EvalReferenceShapeFunctions(refcoords)
   */
  Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override;

  /**
   * @copydoc
   * lf::uscalfe::ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions(refcoords)
   */
  Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override;

  /**
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::EvaluationNodes()
   */
  Eigen::MatrixXd EvaluationNodes() const override;

  /**
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  size_type NumEvaluationNodes() const override;

  /**
   * @copydoc
   * lf::uscalfe::ScalarReferenceFiniteElement::NodalValuesToDofs(nodvals)
   */
  Eigen::Matrix<scalar_type, 1, Eigen::Dynamic> NodalValuesToDofs(
      const Eigen::Matrix<scalar_type, 1, Eigen::Dynamic>& nodvals)
      const override;
};

}  // namespace NonConformingCrouzeixRaviartFiniteElements

#endif  // NUMPDE_CRREFERENCEFINITEELEMENT_H
