#ifndef InitCondLV_CC_
#define InitCondLV_CC_
/**
 * @file initcondlv.cc
 * @brief NPDE homework InitCondLV code
 * @author lfilippo, tille, jgacon, dcasati
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <iostream>
#include <utility>

#include "ode45.h"

namespace InitCondLV {

/* Compute the maps Phi(t,y0) and W(t,y0) at final time T.
 * Use initial data given by u0 and v0. */
/* SAM_LISTING_BEGIN_1 */
std::pair<Eigen::Vector2d, Eigen::Matrix2d> PhiAndW(double u0, double v0,
                                                    double T) {
  // Save the values of Phi and W at time T in PaW.first and PaW.second resp.
  std::pair<Eigen::Vector2d, Eigen::Matrix2d> PaW;

  //====================
  // Your code goes here
  //====================
  return PaW;
}
/* SAM_LISTING_END_1 */

}  // namespace InitCondLV

#endif  // #define InitCondLV_CC_
