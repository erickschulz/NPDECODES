/**
 * @file crossprod.cc
 * @brief NPDE homework CrossProd code
 * @author Unknown, Oliver Rietmann
 * @date 31.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "crossprod.h"

#include <Eigen/Geometry>
#include <iomanip>
#include <iostream>
#include <vector>

namespace CrossProd {

/* SAM_LISTING_BEGIN_0 */
void tab_crossprod() {
  // TO DO (13-1.e): solve the cross-product ODE with the implicit RK method
  // defined in solve_imp_mid. Tabulate the norms of the results at all steps.
#if SOLUTION
  double T = 10.;
  int M = 128;
  // set data
  double c = 1.;
  Eigen::Vector3d y0(1., 1., 1.);
  Eigen::Vector3d a(1., 0., 0.);

  // define rhs
  auto f = [a, c](Eigen::Vector3d y) -> Eigen::Vector3d {
    return a.cross(y) + c * y.cross(a.cross(y));
  };
  // define Jacobian of rhs
  auto Jf = [a, c](Eigen::Vector3d y) -> Eigen::Matrix3d {
    Eigen::Matrix3d temp;
    temp << -c * (a(1) * y(1) + a(2) * y(2)),
        c * (2 * a(0) * y(1) - a(1) * y(0)) - a(2),
        a(1) + c * (2 * a(0) * y(2) - a(2) * y(0)),
        a(2) - c * (a(0) * y(1) - 2 * a(1) * y(0)),
        -c * (a(0) * y(0) + a(2) * y(2)),
        c * (2 * a(1) * y(2) - a(2) * y(1)) - a(0),
        -a(1) - c * (a(0) * y(2) - 2 * a(2) * y(0)),
        a(0) - c * (a(1) * y(2) - 2 * a(2) * y(1)),
        -c * (a(0) * y(0) + a(1) * y(1));
    return temp;
  };

  std::vector<Eigen::VectorXd> res_imp = solve_imp_mid(f, Jf, T, y0, M);

  std::cout << "1. Implicit midpoint method" << std::endl;
  std::cout << std::setw(10) << "t" << std::setw(15) << "norm(y(t))"
            << std::endl;

  for (int i = 0; i < M + 1; ++i) {
    std::cout << std::setw(10) << T * i / M << std::setw(15)
              << res_imp[i].norm() << std::endl;
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  /* SAM_LISTING_END_0 */

  /* SAM_LISTING_BEGIN_1 */
  // TO DO (13-1.g): solve the cross-product ODE with the implicit RK method
  // defined in solve_lin_mid. Tabulate the norms of the results at all steps.
#if SOLUTION
  std::vector<Eigen::VectorXd> res_lin = solve_lin_mid(f, Jf, T, y0, M);
  std::cout << "\n2. Linear implicit midpoint method" << std::endl;
  std::cout << std::setw(10) << "t" << std::setw(15) << "norm(y(t))"
            << std::endl;
  for (int i = 0; i < M + 1; ++i) {
    std::cout << std::setw(10) << T * i / M << std::setw(15)
              << res_lin[i].norm() << std::endl;
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  /* SAM_LISTING_END_1 */
}

}  // namespace CrossProd
