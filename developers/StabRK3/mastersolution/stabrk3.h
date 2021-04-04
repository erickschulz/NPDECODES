#ifndef STABRK3_H_
#define STABRK3_H_

/**
 * @file stabrk3.h
 * @brief NPDE homework StabRK3 code
 * @author Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <vector>

namespace StabRK3 {

Eigen::Vector2d predPrey(Eigen::Vector2d y0, double T, unsigned N);

std::vector<Eigen::Vector2d> simulatePredPrey(
    const std::vector<unsigned int> &N_list);

}  // namespace StabRK3

#endif  // #ifndef STABRK3_H_
