#include <Eigen/Core>
#include <iostream>
#include <utility>

#include "initcondlv.h"

int main() {
  // The test uses the input u0=2.8, v0=1.5, T=2
  std::pair<Eigen::Vector2d, Eigen::Matrix2d> PaW =
      InitCondLV::PhiAndW(2.8, 1.5, 2);
  std::cout << "Test of PhiAndW():\nPhi = " << PaW.first.transpose()
            << "\nW = \n"
            << PaW.second << "\n\n";

  Eigen::Vector2d y = InitCondLV::findInitCond();
  std::cout << "Test of findInitCond():\ny = " << y.transpose() << "\n";
  return 0;
}
