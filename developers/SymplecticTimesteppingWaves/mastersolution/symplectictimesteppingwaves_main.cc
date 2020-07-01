/**
 * @file
 * @brief NPDE homework SymplecticTimesteppingWaves code
 * @author Erick Schulz
 * @date  14/04/2019
 * @copyright Developed at ETH Zurich
 */

#include "symplectictimesteppingwaves.h"
#include "symplectictimesteppingwaves_ode.h"

using namespace SymplecticTimesteppingWaves;

int main(int /*argc*/, char** /*argv*/) {
  // Tabulated study convergence for the symplectic stepping method for the
  // harmonic oscillator ODE dp/dt = q, dq/dt = -p
  sympTimesteppingODETest();

  // Solve the wave equation using the symplectic stepping
  unsigned int m = 2000;
  wavePropSimulation(m);

  double max_step_size = testStab();
  std::cout << "Maximum uniform step size for stability is roughly "
            << max_step_size << std::endl;
}
