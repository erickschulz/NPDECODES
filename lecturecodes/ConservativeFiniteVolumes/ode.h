#ifndef ODE_HPP
#define ODE_HPP

#include <Eigen/Core>

// HACK BEGIN
// Avoid compatibility issue between
// boost/numeric/odeint/external/eigen/eigen_algebra.hpp by shaddwing it with
// modified header ./eigen_algebra.hpp
#include "eigen_algebra.hpp"
// HACK END

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>

#include <tuple>

namespace ConsFV {

/**
 * @brief Tries to mimic ode45 from Matlab
 * Integrates the system of differential equations y' = f(t, y) from time t0 to
 * tfinal with initial conditions y0.
 * @param odefun The function f(t, y) for a scalar t and a vector y.
 * @param t0 The starting time.
 * @param tfinal The ending time.
 * @param y0 The initial conditions.
 * @param abserr The absolute error tolerance.
 * @param relerr The relative error tolerance.
 * @return
 */
inline std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
ode45(const std::function<void(const Eigen::VectorXd &x, Eigen::VectorXd &dxdt,
                               const double t)> &odefun,
      const double t0, const double tfinal, const Eigen::VectorXd &y0,
      const double abserr = 1.0E-8, const double relerr = 1.0E-6) {
  // initialization
  Eigen::VectorXd mu0 = y0;
  Eigen::VectorXd timesteps;
  Eigen::MatrixXd mus(y0.size(), 0);
  // initialize a runge kutta stepper with adaptive step size
  auto stepper = boost::numeric::odeint::make_controlled(
      abserr, relerr,
      boost::numeric::odeint::runge_kutta_dopri5<Eigen::VectorXd>());
  // solve the system recording the time point and corresponding solution vector
  // at every step
  boost::numeric::odeint::integrate_adaptive(
      stepper, odefun, mu0, t0, tfinal, 0.01,
      [&mus, &timesteps](const Eigen::VectorXd &x, double t) {
        mus.conservativeResize(mus.rows(), mus.cols() + 1);
        mus.col(mus.cols() - 1) = x;
        timesteps.conservativeResize(timesteps.size() + 1);
        timesteps(timesteps.size() - 1) = t;
      });
  return std::make_tuple(timesteps, mus);
}

} // namespace ConsFV

#endif
