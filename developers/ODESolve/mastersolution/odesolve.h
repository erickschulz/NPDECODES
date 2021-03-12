#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace ODESolve {

//! \brief Evolves the vector y0 using the evolution operator \tilde{\Psi} with
//! step-size y0 \tparam DiscEvlOp type for evolution operator (e.g. lambda
//! function type) \param[in] Psi original evolution operator, must have
//! operator(double, const Eigen::VectorXd&) \param[in] p parameter p for
//! construction of Psi tilde \param[in] h step-size \param[in] y0 previous step
//! \return Evolved step \f$ \tilde{\Psi}^h y0 \f$
/* SAM_LISTING_BEGIN_0 */
template <class DiscEvlOp>
Eigen::VectorXd psitilde(const DiscEvlOp& Psi, unsigned int p, double h,
                         const Eigen::VectorXd& y0) {
#if SOLUTION
  return (Psi(h, y0) - (2 << (p - 1)) * Psi(h / 2., Psi(h / 2., y0))) /
         (1. - (2 << (p - 1)));
#else   // TEMPLATE
  // TODO: implement psi tilde operator and replace the dummy return value y0
  return y0;
#endif  // TEMPLATE
}
/* SAM_LISTING_END_0 */

//! \brief Evolves the vector y0 using the evolution operator from time 0 to T
//! using equidistant steps \tparam DiscEvlOp type for evolution operator (e.g.
//! lambda function type) \param[in] Psi original evolution operator, must have
//! operator(double, const Eigen::VectorXd&) \param[in] T final time \param[in]
//! y0 initial data \param[in] N number of steps \return Eigen::VectorXd of all
//! steps y_0, y_1, ...
/* SAM_LISTING_BEGIN_1 */
template <class DiscEvlOp>
std::vector<Eigen::VectorXd> odeintequi(const DiscEvlOp& Psi, double T,
                                        const Eigen::VectorXd& y0, int N) {
  double h = T / N;
  double t = 0.;
  std::vector<Eigen::VectorXd> Y;

#if SOLUTION
  Y.reserve(N + 1);
  Y.push_back(y0);
  Eigen::VectorXd y = y0;

  while (t < T) {
    y = Psi(h, Y.back());

    Y.push_back(y);
    t += std::min(T - t, h);
  }
#else   // TEMPLATE
  // TODO: implement equidistant solver for evolution operaotr Psi
#endif  // TEMPLATE

  return Y;
}
/* SAM_LISTING_END_1 */

//! \brief Evolves the vector y0 using the evolution operator from time 0 to T
//! using adaptive error control \tparam DiscEvlOp type for evolution operator
//! (e.g. lambda function type) \param[in] Psi original evolution operator, must
//! have operator(double, const Eigen::VectorXd&) \param[in] T final time
//! \param[in] y0 initial data
//! \param[in] h0 initial step size
//! \param[in] p parameter p for construction of Psi tilde
//! \param[in] reltol relative tolerance for error control
//! \param[in] abstol absolute tolerance for error control
//! \param[in] p parameter p for construction of Psi tilde
//! \param[in] hmin minimal step size
//! \return Eigen::VectorXd of all steps y_0, y_1, ...
/* SAM_LISTING_BEGIN_3 */
template <class DiscEvlOp>
std::vector<Eigen::VectorXd> odeintssctrl(const DiscEvlOp& Psi, unsigned int p,
                                          const Eigen::VectorXd& y0, double T,
                                          double h0, double reltol,
                                          double abstol, double hmin) {
  double t = 0.;
  double h = h0;
  std::vector<Eigen::VectorXd> Y;
  Y.push_back(y0);
  Eigen::VectorXd y = y0;

  while (t < T && h > hmin) {
#if SOLUTION
    Eigen::VectorXd y_high = psitilde(Psi, p, std::min(T - t, h), y);
    Eigen::VectorXd y_low = Psi(std::min(T - t, h), y);
    double est = (y_high - y_low).norm();
    double tol = std::max(reltol * y.norm(), abstol);
    h = h * std::max(0.5, std::min(2., std::pow(tol / est, 1. / (p + 1))));

    if (est < tol) {
#else   // TEMPLATE
    Eigen::VectorXd y_high =
        Eigen::VectorXd::Zero(y0.size());  // TODO: fix this line
    Eigen::VectorXd y_low =
        Eigen::VectorXd::Zero(y0.size());  // TODO: fix this line
    double est = (y_high - y_low).norm();
    // TODO: update $h$

    if (true /* TODO: fix this line */) {
#endif  // TEMPLATE
      y = y_high;
      Y.push_back(y_high);
      t += std::min(T - t, h);
    }
  }
  if (h < hmin) {
    std::cerr << "Warning: Failure at t=" << t
              << ". Unable to meet integration tolerances without reducing the "
                 "step size below the smallest value allowed ("
              << hmin << ") at time t." << std::endl;
  }

  return Y;
}
/* SAM_LISTING_END_3 */

}  // namespace ODESolve
