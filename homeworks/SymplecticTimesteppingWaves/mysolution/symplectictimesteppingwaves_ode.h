#ifndef SYMPLECTIC_ODE_HPP
#define SYMPLECTIC_ODE_HPP

/**
 * @file
 * @brief NPDE homework SymplecticTimesteppingWaves
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

#include <cmath>
#include <iostream>

#include <Eigen/Dense>

namespace SymplecticTimesteppingWaves {

void sympTimestep(double tau, Eigen::Vector2d &pq_j);
Eigen::Vector2d sympTimesteppingHarmonicOscillatorODE(unsigned int m);
void sympTimesteppingODETest();

} // namespace SymplecticTimesteppingWaves
#endif
