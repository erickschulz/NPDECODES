/**
 * @file stabrk3.cc
 * @brief NPDE homework StabRK3 code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "stabrk3.h"

#include <Eigen/Core>
#include <vector>

#include "rkintegrator.h"

namespace StabRK3 {

/* SAM_LISTING_BEGIN_0 */
Eigen::Vector2d predPrey(Eigen::Vector2d y0, double T, unsigned int N) {
  double h = T / N;
  Eigen::Vector2d y = y0;

  //====================
  // Your code goes here
  //====================

  return y;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::Vector2d> simulatePredPrey(
    const std::vector<unsigned int> &N_list) {
  int M = N_list.size();
  std::vector<Eigen::Vector2d> yT_list(M);

  //====================
  // Your code goes here
  //====================

  return yT_list;
}
/* SAM_LISTING_END_1 */

}  // namespace StabRK3
