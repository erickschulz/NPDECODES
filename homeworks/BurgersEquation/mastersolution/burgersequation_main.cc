/**
 * @file burgersequation_main.cc
 * @brief NPDE homework BurgersEquation code
 * @author Oliver Rietmann
 * @date 15.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "burgersequation.h"

#include <fstream>
#include <iostream>

#include <Eigen/Core>

int main() {
  /* SAM_LISTING_BEGIN_1 */
  const unsigned int N = 100;
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N + 1, -1.0, 4.0);
  Eigen::VectorXd mu03 = BurgersEquation::solveBurgersGodunov(0.3, N);
  Eigen::VectorXd mu30 = BurgersEquation::solveBurgersGodunov(3.0, N);

  // Write the solutions to a file that can be used for plotting.
  /* SOLUTION_BEGIN */
  std::ofstream solution_file;
  solution_file.open("solution.csv");
  solution_file << x.transpose().format(BurgersEquation::CSVFormat)
                << std::endl;
  solution_file << mu03.transpose().format(BurgersEquation::CSVFormat)
                << std::endl;
  solution_file << mu30.transpose().format(BurgersEquation::CSVFormat)
                << std::endl;
  solution_file.close();
  std::cout
      << "Generated solution.csv; now run plot_solution_mastersolution.py."
      << std::endl;
  /* SOLUTION_END */
  /* SAM_LISTING_END_1 */

  /* SAM_LISTING_BEGIN_2 */
  Eigen::Matrix<double, 3, 4> result = BurgersEquation::numexpBurgersGodunov();

  // Write the result to a file that can be used for plotting.
  /* SOLUTION_BEGIN */
  std::ofstream error_file;
  error_file.open("error.csv");
  error_file << result.format(BurgersEquation::CSVFormat) << std::endl;
  error_file.close();
  std::cout << "Generated error.csv; now run plot_error_mastersolution.py."
            << std::endl;
  /* SOLUTION_END */
  /* SAM_LISTING_END_2 */

  return 0;
}
