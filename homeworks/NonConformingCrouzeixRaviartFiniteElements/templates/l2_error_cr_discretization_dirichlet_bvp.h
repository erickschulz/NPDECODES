/**
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss, edited Amélie Loher
 * @date   18.03.2019, 03.03.20
 * @copyright Developed at ETH Zurich
 */

#ifndef NUMPDE_L2_ERROR_CR_DISCRETIZATION_DIRICHLET_BVP_H
#define NUMPDE_L2_ERROR_CR_DISCRETIZATION_DIRICHLET_BVP_H

#include <cmath>
#include <string>

#include <lf/io/io.h>
#include <lf/mesh/mesh.h>

#include "cr_fe_space.h"
#include "cr_l2_error.h"
#include "solve_cr_dirichlet_bvp.h"

namespace NonConformingCrouzeixRaviartFiniteElements {

double L2errorCRDiscretizationDirichletBVP(const std::string &filename) {
  double l2_error;

// TODO: task 2-14.x)
  //====================
  // Your code goes here
  //====================
  return l2_error;
}

}  // namespace NonConformingCrouzeixRaviartFiniteElements

#endif  // NUMPDE_L2_ERROR_CR_DISCRETIZATION_DIRICHLET_BVP_H
