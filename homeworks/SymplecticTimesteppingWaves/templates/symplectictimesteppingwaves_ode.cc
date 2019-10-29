/**
 * @file
 * @brief NPDE homework SymplecticTimesteppingWaves
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "symplectictimesteppingwaves_ode.h"

namespace SymplecticTimesteppingWaves {

void sympTimestep(double tau, Eigen::Vector2d& pq_j) {
  /* SOLUTION_BEGIN */
  /* TODO Your implementation goes here! */
  /* SOLUTION_END */
}

Eigen::Vector2d sympTimesteppingHarmonicOscillatorODE(unsigned int m) {
  Eigen::Vector2d approx_sol;
  approx_sol << 0.0, 1.0;  // initial conditions
  double tau = 2.0 * M_PI / m;
  for (int i = 0; i < m; i++) {
    sympTimestep(tau, approx_sol);
  }
  return approx_sol;
}

void sympTimesteppingODETest() {
  /* SOLUTION_BEGIN */
  /* TODO Your implementation goes here! */
  /* SOLUTION_END */
}

}  // namespace SymplecticTimesteppingWaves
