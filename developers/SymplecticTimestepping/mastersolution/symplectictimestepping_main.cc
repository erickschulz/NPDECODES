#include <Eigen/Dense>
#include <iostream>

#include "symplectictimestepping.h"

int main() {
  Eigen::Vector2d pq0, pq1;
  pq0 << 0., 1.;
  pq1 << 0.96357170038053707728, 0.12823014413913069731;
  double tol = 1E-8;
  SymplecticTimestepping::sympTimestep(0.5 * SymplecticTimestepping::PI, pq0);

  std::cout << "Test for 'sympTimestep()': ";
  if ((pq0 - pq1).norm() < tol) {
    std::cout << "passed!\n";
  } else {
    std::cout << "failed\n";
  }

  SymplecticTimestepping::sympTimesteppingODETest();
  std::cout << "No test for 'sympTimesteppingODETest()'.\n";

  Eigen::Vector3d p0, q0;
  p0 << 0.1, 0.0, 0.5;
  q0 << 0.0, 0.1, 0.5;
  int M = 5;
  double T = 1.1;
  Eigen::MatrixXd PQsol(3 * 2, M + 1);
  PQsol << 0.1, 0.09678001159252096, 0.08529466710449594, 0.06584032145575656,
      0.04296862323725763, 0.02256010635860073, 0, -0.02738875930643383,
      -0.0602025556404987, -0.09072015056737635, -0.1119069922852328,
      -0.1218541432646306, 0.5, 0.3469562614304357, 0.1254605573199863,
      -0.1243991455580988, -0.3446918452398756, -0.4964701845301493, 0,
      0.02176457797371695, 0.04196541799960365, 0.05873659400349878,
      0.07074757360641674, 0.07788776841519347, 0.1, 0.09716774215801476,
      0.08762065486164657, 0.0709505549927613, 0.04847378550452989,
      0.02256428674940955, 0.5, 0.5946616006586585, 0.6479303643062508,
      0.6484357449813002, 0.596106795554733, 0.502260275823015;
  tol = 1E-4;
  std::cout << "Test for 'simulateHamiltonianDynamics()': ";
  if ((PQsol -
       SymplecticTimestepping::simulateHamiltonianDynamics(p0, q0, T, M))
          .norm() < tol) {
    std::cout << "passed!" << std::endl;
  } else {
    std::cout << "failed!" << std::endl;
  }

  return 0;
}
