#include <Eigen/Core>
#include <Eigen/QR>
#include <iomanip>
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>

#include "ode45.h"

namespace MatODE {

//! \brief Solve matrix IVP Y' = -(Y-Y')*Y using ode45 up to time T
//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
//! \return Matrix of solution of IVP at t = T
/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd matode(const Eigen::MatrixXd& Y0, double T) {
  auto F = [](const Eigen::MatrixXd& M) { return -(M - M.transpose()) * M; };
  ode45<Eigen::MatrixXd> O(F);

  // Set tolerances
  O.options.atol = 10e-10;
  O.options.rtol = 10e-8;

  // Return only matrix at $T$, (solution is vector
  // of pairs $(y(t_k), t_k)$ for each step k
  return O.solve(Y0, T).back().first;
}
/* SAM_LISTING_END_1 */

//! \brief Find if invariant is preserved after evolution with matode
//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
//! \return true if invariant was preserved (up to round-off),
//! i.e. if norm was less than 10*eps
/* SAM_LISTING_BEGIN_2 */
bool checkinvariant(const Eigen::MatrixXd& M, double T) {
  Eigen::MatrixXd N(3, 3);

  N = matode(M, T);

  if ((N.transpose() * N - M.transpose() * M).norm() <
      10 * std::numeric_limits<double>::epsilon() * M.norm()) {
    return true;
  } else {
    return false;
  }
}
/* SAM_LISTING_END_2 */

//! \brief Implement ONE step of explicit Euler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd expeulstep(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Y0,
                           double h) {
#if SOLUTION
  return Y0 + h * A * Y0;
#else   // TEMPLATE
  // TODO: implement one explicit Euler step
  return Y0;
#endif  // TEMPLATE
}
/* SAM_LISTING_END_3 */

//! \brief Implement ONE step of implicit Euler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
/* SAM_LISTING_BEGIN_4 */
Eigen::MatrixXd impeulstep(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Y0,
                           double h) {
#if SOLUTION
  return (Eigen::MatrixXd::Identity(3, 3) - h * A).partialPivLu().solve(Y0);
#else   // TEMPLATE
  // TODO: implement one implicit Euler step
  return Y0;
#endif  // TEMPLATE
}
/* SAM_LISTING_END_4 */

//! \brief Implement ONE step of implicit midpoint ruler applied to Y0, of ODE
//! Y' = A*Y \param[in] A matrix A of the ODE \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
/* SAM_LISTING_BEGIN_5 */
Eigen::MatrixXd impstep(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Y0,
                        double h) {
#if SOLUTION
  return (Eigen::MatrixXd::Identity(3, 3) - h * 0.5 * A)
      .partialPivLu()
      .solve(Y0 + h * 0.5 * A * Y0);
#else   // TEMPLATE
  // TODO: implement one implicit midpoint step.
  return Y0;
#endif  // TEMPLATE
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
std::tuple<double, double, double> checkOrthogonality(void) {
  // TO DO (12-1.c): compute and tabulate the Frobenius norms of Y_k'*Y_k - I
  // for 20 steps of eeulstep, ieulstep and impstep.
  // Return the values corresponding to the last step.
  unsigned int n = 3;
  double h = 0.01;
  Eigen::MatrixXd M(n, n);
  M << 8, 1, 6, 3, 5, 7, 9, 9, 2;
  std::vector<double> norms(3);

#if SOLUTION
  // Build Q
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(n, n);
  qr.compute(M);
  Eigen::MatrixXd Q = qr.householderQ();

  // Build A
  Eigen::MatrixXd A(n, n);
  A << 0, 1, 1, -1, 0, 1, -1, -1, 0;
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);

  Eigen::MatrixXd Meeul = Q, Mieul = Q, Mimp = Q;

  std::vector<int> sep = {5, 15};
  std::cout << std::setw(sep[0]) << "step" << std::setw(sep[1]) << "exp. Eul"
            << std::setw(sep[1]) << "imp. Eul" << std::setw(sep[1]) << "Mid-Pt"
            << std::endl;
  // Norm of Y'Y-I for initial value
  std::cout << std::setw(sep[0]) << "-1" << std::setw(sep[1])
            << (Meeul.transpose() * Meeul - I).norm() << std::setw(sep[1])
            << (Mieul.transpose() * Mieul - I).norm() << std::setw(sep[1])
            << (Mimp.transpose() * Mimp - I).norm() << std::endl;

  // Norm of Y'Y-I for 20 steps
  for (unsigned int j = 0; j < 20; ++j) {
    Meeul = expeulstep(A, Meeul, h);
    Mieul = impeulstep(A, Mieul, h);
    Mimp = impstep(A, Mimp, h);

    norms[0] = (Meeul.transpose() * Meeul - I).norm();
    norms[1] = (Mieul.transpose() * Mieul - I).norm();
    norms[2] = (Mimp.transpose() * Mimp - I).norm();

    std::cout << std::setw(sep[0]) << j << std::setw(sep[1]) << norms[0]
              << std::setw(sep[1]) << norms[1] << std::setw(sep[1]) << norms[2]
              << std::endl;
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

  return std::make_tuple(norms[0], norms[1], norms[2]);
}
/* SAM_LISTING_END_6 */

}  // namespace MatODE
