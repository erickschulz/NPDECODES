/**
 * @ file
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author Tobias Rohner
 * @ date 25-03-2022
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "mehrstellenverfahren.h"

#include <iomanip>
#include <iostream>
#include <vector>

namespace mehrstellenverfahren {

/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<double> compMehrstellenA(unsigned int M) {
  // For the sake of efficiency the use of Eigen's sparse matrix data type is
  // essential. The matrix is stored in CCS format.
  Eigen::SparseMatrix<double> A(M * M, M * M);
  //====================
  // Your code goes here
  //====================
  return A;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double compgriderr(unsigned int M) {
  double err = 0;
  //====================
  // Your code goes here
  //====================
  return err;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void tabulateMehrstellenError() {
  std::vector<double> errs;
  std::vector<unsigned int> Ms = {5, 10, 20, 40, 80, 160};
  //====================
  // Your code goes here
  //====================
  // Print the collected data out in a table
  std::cout << "| M   | h       | error        |\n";
  std::cout << "-------------------------------|\n";
  for (int i = 0; i < Ms.size(); ++i) {
    const unsigned int M = Ms[i];
    const double h = 1. / (M + 1);
    const double err = errs[i];
    std::cout << "| " << std::setw(3) << M << " | " << std::fixed
              << std::setprecision(5) << h << " | " << std::scientific << err
              << " |\n";
  }
}
/* SAM_LISTING_END_3 */

}  // namespace mehrstellenverfahren
