/**
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss
 * @date   18.03.2019
 * @copyright Developed at ETH Zurich
 */

#ifndef NUMPDE_SOLVE_CR_DIRICHLET_BVP_H
#define NUMPDE_SOLVE_CR_DIRICHLET_BVP_H

#include <lf/assemble/assemble.h>
#include <lf/uscalfe/uscalfe.h>
#include "cr_fe_space.h"

namespace NonConformingCrouzeixRaviartFiniteElements {

/* SAM_LISTING_BEGIN_1 */
template <typename GAMMA_COEFF, typename F_FUNCTOR>
Eigen::VectorXd solveCRDirichletBVP(std::shared_ptr<CRFeSpace> fe_space,
                                    GAMMA_COEFF &&gamma, F_FUNCTOR &&f) {
  Eigen::VectorXd sol;

  // TODO: task 2-14.v)
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */

  return sol;
}
/* SAM_LISTING_END_1 */

}  // namespace NonConformingCrouzeixRaviartFiniteElements

#endif  // NUMPDE_SOLVE_CR_DIRICHLET_BVP_H
