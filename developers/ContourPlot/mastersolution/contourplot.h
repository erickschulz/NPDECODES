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

  // Right-hand-side vector field of isoline ODE
  auto rhs = [gradF](Eigen::Vector2d x) -> Eigen::Vector2d {
    Eigen::Vector2d gradFx = gradF(x);
    return Eigen::Vector2d(-gradFx(1), gradFx(0)) / gradFx.norm();
  };

#if SOLUTION
  // Adaptive explicit Runge-Kutta method
  ode45<Eigen::Vector2d> integrator(rhs);
  // Set tolerances for timestep control
  integrator.options.rtol = 0.00001;
  integrator.options.atol = 1e-9;
  // Perform explicit timestepping
  std::vector<std::pair<Eigen::Vector2d, double>> sol = integrator.solve(y0, T);
  // Convert output into requested format: points on isoline arranged into the
  // columns of a matrix.
  int M = sol.size() - 1;
  States = Eigen::MatrixXd::Zero(2, M + 1);
  for (int m = 0; m <= M; ++m) {
    States.col(m) = sol[m].first;
  }
#else
  //====================
  // Your code goes here
  // Use ode45 to solve the ODE with right-hand side
  // given by the function rhs(...) and fill the
  // matrix States with the solution vectors
  //====================
#endif

  return States;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 2, Eigen::Dynamic> crookedEgg() {
#if SOLUTION
  auto gradF = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return 4.0 * x.squaredNorm() * x - 3.0 * x.cwiseAbs2();
  };
  Eigen::Vector2d y0(1.0, 0.0);
  double T = 4.0;
  return computeIsolinePoints(gradF, y0, T);
#else
  //====================
  // Your code goes here
  // Replace the following dummy return value
  // by the matrix containing the isoline points:
  return Eigen::Matrix<double, 2, 42>::Zero();
  //====================
#endif
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <typename FFunctor>
Eigen::Matrix<double, 2, Eigen::Dynamic> computeIsolinePointsDQ(
    FFunctor &&F, Eigen::Vector2d y0, double T) {

#if SOLUTION
  // Get gradient of F by symmetric difference quotients with
  // a span given by the root of the machine precision in order
  // to curb cancellation as much as possible.
  double h = std::sqrt(std::numeric_limits<double>::epsilon());
  auto gradF = [h, F](Eigen::Vector2d x) -> Eigen::Vector2d {
    Eigen::Vector2d dx(h, 0.0);
    Eigen::Vector2d dy(0.0, h);
    return Eigen::Vector2d(F(x + dx) - F(x - dx), F(x + dy) - F(x - dy)) / (2.0 * h);
  };

  return computeIsolinePoints(gradF, y0, T);
#else
  //====================
  // Your code goes here
  // Replace the following dummy return value
  // by the matrix containing the isoline points:
  return Eigen::Matrix<double, 2, 42>::Zero();
  //====================
#endif
}
/* SAM_LISTING_END_2 */

}  // namespace ContourPlot
