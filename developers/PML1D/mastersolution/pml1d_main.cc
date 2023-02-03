/**
 * @ file pml1d_main.cc
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author R. Hiptmair
 * @ date January 2023
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <iostream>

#include "pml1d.h"

int main(int /*argc*/, char ** /*argv*/) {
  PML1D::tabulateExp1();
  PML1D::plotExp(200, 200, 4.0);
  return 0;
}
