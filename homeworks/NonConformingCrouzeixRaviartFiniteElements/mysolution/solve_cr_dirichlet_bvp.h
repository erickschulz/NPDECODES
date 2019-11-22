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

namespace NonConformingCrouzeixRaviartFiniteElements
{

template <typename GAMMA_COEFF, typename F_FUNCTOR>
Eigen::VectorXd solveCRDirichletBVP(std::shared_ptr<CRFeSpace> fe_space,
                                    GAMMA_COEFF &&gamma, F_FUNCTOR &&f)
{
    Eigen::VectorXd sol;
  //====================
  // Your code goes here
  //====================
    return sol;
}

} // namespace NonConformingCrouzeixRaviartFiniteElements

#endif // NUMPDE_SOLVE_CR_DIRICHLET_BVP_H
