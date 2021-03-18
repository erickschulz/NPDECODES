#ifndef LV_H_
#define LV_H_
/**
 * @file initcondlv.cc
 * @brief NPDE homework InitCondLV code
 * @author lfilippo, tille, jgacon, dcasati
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

namespace InitCondLV {

std::pair<Eigen::Vector2d, Eigen::Matrix2d> PhiAndW(double u0, double v0,
                                                    double T);

}  // namespace InitCondLV

#endif  // #define LV_H_
