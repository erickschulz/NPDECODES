#ifndef SYMPLECTIC_ODE_HPP
#define SYMPLECTIC_ODE_HPP

/**
 * @file
 * @brief NPDE homework SymplecticTimesteppingWaves
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include <Eigen/Dense>
#include <iostream>

namespace SymplecticTimesteppingWaves {

void sympTimestep(double tau, Eigen::Vector2d& pq_j);
Eigen::Vector2d sympTimesteppingHarmonicOscillatorODE(unsigned int m);
void sympTimesteppingODETest();

}  // namespace SymplecticTimesteppingWaves
#endif
