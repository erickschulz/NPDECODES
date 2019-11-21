/** @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss
 * @date   16.03.2019
 * @copyright Developed at ETH Zurich */

#ifndef CRREFERENCEFINITEELEMENT_H
#define CRREFERENCEFINITEELEMENT_H

#include <lf/base/base.h>
#include <lf/uscalfe/uscalfe.h>
#include "cr_types.h"

namespace NonConformingCrouzeixRaviartFiniteElements {

class CRReferenceFiniteElement final
    : public lf::uscalfe::ScalarReferenceFiniteElement<scalar_type> {
 public:
  lf::base::RefEl RefEl() const override;
  unsigned int Degree() const override;
  size_type NumRefShapeFunctions() const override;
  size_type NumRefShapeFunctions(dim_t codim) const override;
  size_type NumRefShapeFunctions(dim_t codim, sub_idx_t subidx) const override;
  Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override;
  Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override;
  Eigen::MatrixXd EvaluationNodes() const override;
  size_type NumEvaluationNodes() const override;
  Eigen::Matrix<scalar_type, 1, Eigen::Dynamic> NodalValuesToDofs(
      const Eigen::Matrix<scalar_type, 1, Eigen::Dynamic>& nodvals)
      const override;
};

}  // namespace NonConformingCrouzeixRaviartFiniteElements

#endif  // CRREFERENCEFINITEELEMENT_H
