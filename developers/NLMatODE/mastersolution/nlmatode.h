#ifndef NLMATODE_H_
#define NLMATODE_H_

#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/LU>

#include "ode45.h"
#include "polyfit.h"


namespace NLMatODE {

//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd matode(const Eigen::MatrixXd &Y0, double T) {
  // TO DO (12-5.a): use the ode45 class to find an approximation
  // of the matrix IVP $Y' = -(Y-Y')*Y$ at time $T$
  Eigen::MatrixXd YT;
  // START
  // Define the RHS
  auto F = [] (const Eigen::MatrixXd & M) {
      return -(M  - M.transpose())*M;
  };
  ode45<Eigen::MatrixXd> O(F);

  // Set tolerances
  O.options.atol = 10e-10;
  O.options.rtol = 10e-8;

  // Return only matrix at $T$, (solution is vector
  // of pairs $(y(t_k), t_k)$ for each step k
  YT = O.solve(Y0, T).back().first;
  // END
  return YT;
}
/* SAM_LISTING_END_1 */


//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
/* SAM_LISTING_BEGIN_2 */
bool checkinvariant(const Eigen::MatrixXd &M, double T) {
  // TO DO (12-5.c): check if $Y'*Y$ is preserved at the time $T$ by matode.
  // START
  Eigen::MatrixXd N = matode(M, T);

  if( (N.transpose()*N-M.transpose()*M).norm() <
      10 * std::numeric_limits<double>::epsilon()* M.norm()) {
      return true;
  } else {
      return false;
  }
  // END
}
/* SAM_LISTING_END_2 */


/* SAM_LISTING_BEGIN_3 */
double cvgDiscreteGradientMethod(void) {
  // TO DO (12-5.d): compute the fitted convergence rate of the Discrete gradient
  // method. Also tabulate the values M and the errors. 
  double conv_rate;
  // START
  double T = 1.0;
  // initial value
  Eigen::MatrixXd Y0 = Eigen::MatrixXd::Zero(5,5);
  Y0(4,0) = 1;
  for (unsigned int i = 0; i < 4; ++i) { 
    Y0(i,i + 1) = 1; 
  }
  // reference solution
  Eigen::MatrixXd Y_ex = matode(Y0, T);
  
  // define the rhs
  auto F = [] (const Eigen::MatrixXd & M) {
      return -(M  - M.transpose())*M;
  };
  
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(5,5);
  Eigen::ArrayXd MM(8);
  Eigen::ArrayXd err(8);
  
  std::cout << "Error for equidistant steps:" << std::endl;
  std::cout << "M" << "\t" << "Error" << std::endl;
  for (unsigned int i = 0; i < 8; ++i) {
    unsigned int M = 10*std::pow(2,i);
    double h = T/M;
    MM(i) = M;
    Eigen::MatrixXd Y = Y0;
    for (unsigned int j = 0; j < M; ++j ) {
      Eigen::MatrixXd Ystar = Y + 0.5 * h * F(Y);
      
      Eigen::MatrixXd Yinc = 0.5 * h * (Ystar - Ystar.transpose());
      Y = (I + Yinc).lu().solve((I - Yinc)*Y);
    }
    err(i) = (Y - Y_ex).norm();
    std::cout << M << "\t" << err(i) << std::endl;
  }
  
  // compute fitted rate
  Eigen::VectorXd coeffs = polyfit(MM.log(),err.log(),1);
  conv_rate = -coeffs(0);
  // END
  return conv_rate;
}
/* SAM_LISTING_END_3 */

}  // namespace NLMatODE

#endif  // #define NLMATODE_H_
