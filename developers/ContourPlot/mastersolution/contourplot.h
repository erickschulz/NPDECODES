#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "ode45.h"

namespace ContourPlot {

/* SAM_LISTING_BEGIN_0 */
template <typename GradientFunctor>
Eigen::Matrix<double, 2, Eigen::Dynamic> computeIsolinePoints(
    GradientFunctor &&gradF, Eigen::Vector2d y0, double T) {
  Eigen::Matrix<double, 2, Eigen::Dynamic> States;
  // To do: (0-2.b)
  // START
  // Right-hand-side vector field of isoline ODE
  auto rhs = [gradF](const Eigen::Vector2d &x) {
    Eigen::Vector2d gradFx = gradF(x);
    Eigen::Vector2d out(-gradFx(1), gradFx(0));  // rotate
    out /= gradFx.norm();                        // scale
    return out;
  };
  // Adaptive explicit Runge-Kutta method
  ode45<Eigen::Vector2d> integrator(rhs);
  // Set tolerances for timestep control
  integrator.options.rtol = 0.00001;
  integrator.options.atol = 1e-9;
  // Perform explicit timestepping
  std::vector<std::pair<Eigen::Vector2d, double>> sol = integrator.solve(y0, T);
  // Convert output into requested format: points on isoline arranged into the
  // columns of a matrix.
  int N = sol.size() - 1;
  States = Eigen::MatrixXd::Zero(2, N + 1);
  for (int i = 0; i <= N; i++) {
    States.col(i) = sol[i].first;
  }
  // END
  return States;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 2, Eigen::Dynamic> crookedEgg() {
  // crookedEggCurve will need to be reshaped to 2*(N+1).
  Eigen::Matrix<double, 2, Eigen::Dynamic> crookedEggCurve;
  // To do: (0-2.c)
  // START
  Eigen::Vector2d y0(1., 0.);
  double T = 4.;
  auto gradF = [](const Eigen::Vector2d &x) {
    double x2 = x(0) * x(0), y2 = x(1) * x(1);
    double x2y2 = 4. * (x2 + y2);
    Eigen::Vector2d grad(x2y2 * x(0) - 3. * x2, x2y2 * x(1) - 3. * y2);
    return grad;
  };
  crookedEggCurve = computeIsolinePoints(gradF, y0, T);
  // END
  return crookedEggCurve;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <typename FFunctor>
Eigen::Matrix<double, 2, Eigen::Dynamic> computeIsolinePointsDQ(
    FFunctor &&F, Eigen::Vector2d y0, double T) {
  // States will need to be reshaped to 2*(N+1).
  Eigen::Matrix<double, 2, Eigen::Dynamic> States;
  // To do: (0-2.d)
  // START
  // Get gradient of F by symmetric difference quotients with
  // a span given by the root of the machine precision in order
  // to curb cancellation as much as possible.
  double epsrt = std::sqrt(std::numeric_limits<double>::epsilon());
  auto gradF = [epsrt, F](const Eigen::Vector2d &x) {
    Eigen::Vector2d grad, dx, dy;
    dx << epsrt, 0.;
    dy << 0., epsrt;
    grad(0) = (F(x + dx) - F(x - dx)) / (2 * epsrt);
    grad(1) = (F(x + dy) - F(x - dy)) / (2 * epsrt);
    return grad;
  };
  // The remainder as in computeIsolinePoints()
  auto rhs = [gradF](const Eigen::Vector2d &x) {
    Eigen::Vector2d gradFx = gradF(x);
    Eigen::Vector2d out(-gradFx(1), gradFx(0));  // rotate
    out /= gradFx.norm();                        // scale
    return out;
  };

  ode45<Eigen::Vector2d> integrator(rhs);
  integrator.options.rtol = 0.00001;
  integrator.options.atol = 1e-9;

  std::vector<std::pair<Eigen::Vector2d, double>> sol = integrator.solve(y0, T);
  int N = sol.size() - 1;
  States = Eigen::MatrixXd::Zero(2, N + 1);
  for (int i = 0; i <= N; i++) {
    States.col(i) = sol[i].first;
  }
  // END
  return States;
}
/* SAM_LISTING_END_2 */

}  // namespace ContourPlot
