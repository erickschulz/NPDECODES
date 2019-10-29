/**
 * @file qfe_provider_tester.h
 * @brief NPDE homework DebuggingFEM code
 * @author Oliver Rietmann
 * @date 03/04/2019
 * @copyright Developed at ETH Zurich
 */

#ifndef QFE_PROV
#define QFE_PROV

#include "qfe_interpolator.h"

#include <lf/assemble/assemble.h>

namespace DebuggingFEM {

template <typename ENTITY_MATRIX_PROVIDER>
class QFEProviderTester {
 public:
  /**
   * @brief Sets up and stores the Galerkin matrix
   * @param dofh dof handler
   * @param element_matrix_provider an instance of LocalLaplaceQFEX, X=1,2,3
   */
  QFEProviderTester(lf::assemble::DofHandler &dofh,
                    ENTITY_MATRIX_PROVIDER &element_matrix_provider);

  /**
   * @brief Computes the energy (H_1-seminorm) of u
   * @param u function of type double(Eigen::Vector2d)
   */
  template <typename FUNCTOR>
  double energyOfInterpolant(FUNCTOR &&u) const;

 private:
  lf::assemble::DofHandler &dofh_;
  ENTITY_MATRIX_PROVIDER &element_matrix_provider_;
  /* Further private data members */
};

template <typename ENTITY_MATRIX_PROVIDER>
QFEProviderTester<ENTITY_MATRIX_PROVIDER>::QFEProviderTester(
    lf::assemble::DofHandler &dofh,
    ENTITY_MATRIX_PROVIDER &element_matrix_provider)
    : dofh_(dofh), element_matrix_provider_(element_matrix_provider) {
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
}

template <typename ENTITY_MATRIX_PROVIDER>
template <typename FUNCTOR>
double QFEProviderTester<ENTITY_MATRIX_PROVIDER>::energyOfInterpolant(
    FUNCTOR &&u) const {
  double energy;
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return energy;
}

}  // namespace DebuggingFEM

#endif
