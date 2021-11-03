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

#include "../../../lecturecodes/Ode45/ode45.h"

namespace InitCondLV {

/* Compute the maps Phi(t,y0) and W(t,y0) at final time T.
 * Use initial data given by u0 and v0. */
/* SAM_LISTING_BEGIN_1 */
std::pair<Eigen::Vector2d, Eigen::Matrix2d> PhiAndW(double u0, double v0,
                                                    double T) {
  // Save the values of Phi and W at time T in PaW.first and PaW.second resp.
  std::pair<Eigen::Vector2d, Eigen::Matrix2d> PaW;

  auto f = [](const Eigen::VectorXd& w) {
    Eigen::VectorXd temp(6);
    temp(0) = (2. - w(1)) * w(0);
    temp(1) = (w(0) - 1.) * w(1);
    temp(2) = (2. - w(1)) * w(2) - w(0) * w(3);
    temp(3) = w(1) * w(2) + (w(0) - 1.) * w(3);
    temp(4) = (2. - w(1)) * w(4) - w(0) * w(5);
    temp(5) = w(1) * w(4) + (w(0) - 1.) * w(5);
    return temp;
  };

  Eigen::VectorXd w0(6);
  w0 << u0, v0, 1., 0, 0, 1.;

  // Construct ode solver with r.h.s
  Ode45<Eigen::VectorXd> O(f);
  // Set options
  O.options.rtol = 1e-14;
  O.options.atol = 1e-12;
  // Solve ODE
  auto sol = O.solve(w0, T);
  // Extract needed component
  Eigen::VectorXd wT = sol.back().first;

  PaW.first << wT(0), wT(1);
  PaW.second << wT(2), wT(4), wT(3), wT(5);
  return PaW;
}
/* SAM_LISTING_END_1 */

}  // namespace InitCondLV

#endif  // #define InitCondLV_CC_
