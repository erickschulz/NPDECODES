#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "contourplot.h"

int main() {
  double tol = 1e-4;
  // Function for crooked egg:
  auto F = [](const Eigen::Vector2d &x) {
    double x2 = x(0) * x(0), y2 = x(1) * x(1);
    double x2y2 = x2 + y2;
    return x2y2 * x2y2 - x2 * x(0) - y2 * x(1);
  };
  // Test function:
  auto G = [](const Eigen::Vector2d &x) {
    return 0.5 * x(0) * x(0) * x(1) - x(1) * x(1) * x(1) / 6;
  };
  // Gradient of test function:
  auto gradG = [](const Eigen::Vector2d &x) {
    Eigen::Vector2d grad(x(0) * x(1), 0.5 * x(0) * x(0) - 0.5 * x(1) * x(1));
    return grad;
  };
  Eigen::Vector2d y0(-0.5, 0.5);
  Eigen::Vector2d z0(1.0, 0.0);

  /*
   *  run computeIsolinePoints()
   */
  Eigen::MatrixXd isoline = ContourPlot::computeIsolinePoints(gradG, y0, 5.);
  std::cout << "Output of computeIsolinePoints() is a " << isoline.rows()
            << " by " << isoline.cols() << " matrix";
  if (isoline.cols() > 5)
    std::cout << ", starting with\n" << isoline.leftCols(5) << " ...\n";
  else
    std::cout << ":\n" << isoline << std::endl;

  /*
   *  test ContourPlot::computeIsolinePoints()
   */
  if (isoline.rows() == 2 && isoline.cols() >= 10) {
    // ode45 takes at least 10 steps by default.
    // Check initial state:
    bool passed_computeIsolinePoints = ((y0 - isoline.col(0)).norm() < tol);
    double Gy0 = G(y0);
    for (int i = 1; i < isoline.cols(); i++) {
      // Check that point lies on the G(y0)-isoline:
      double err = std::abs(G(isoline.col(i)) - Gy0);
      // Check that point is not stationary:
      double step = (isoline.col(i) - isoline.col(i - 1)).norm();
      if (step < tol * tol || err > tol) {
        passed_computeIsolinePoints = false;
        break;
      }
    }
    if (passed_computeIsolinePoints)
      std::cout << "Test for computeIsolinePoints() passed!\n\n";
    else
      std::cout << "Test for computeIsolinePoints() failed.\n\n";
  } else
    std::cout << "Test for computeIsolinePoints() failed: wrong size.\n\n";

  /*
   *  run ContourPlot::crookedEgg()
   */
  Eigen::MatrixXd crookedEggCurve = ContourPlot::crookedEgg();
  std::cout << "Output of crookedEgg() is a " << crookedEggCurve.rows()
            << " by " << crookedEggCurve.cols() << " matrix";
  if (crookedEggCurve.cols() > 5)
    std::cout << ", starting with\n" << crookedEggCurve.leftCols(5) << " ...\n";
  else
    std::cout << ":\n" << crookedEggCurve << std::endl;

  /*
   *  test ContourPlot::crookedEgg()
   */
  if (crookedEggCurve.rows() == 2 && crookedEggCurve.cols() >= 10) {
    // ode45 takes at least 10 steps by default.
    // Check initial state:
    bool passed_crookedEgg = ((z0 - crookedEggCurve.col(0)).norm() < tol);
    for (int i = 1; i < crookedEggCurve.cols(); i++) {
      // Check that point lies on the 0-isoline:
      double err = std::abs(F(crookedEggCurve.col(i)));
      // Check that point is not stationary:
      double step =
          (crookedEggCurve.col(i) - crookedEggCurve.col(i - 1)).norm();
      if (step < tol * tol || err > tol) {
        passed_crookedEgg = false;
        break;
      }
    }
    if (passed_crookedEgg)
      std::cout << "Test for crookedEgg() passed!\n\n";
    else
      std::cout << "Test for crookedEgg() failed.\n\n";
  } else
    std::cout << "Test for crookedEgg() failed: wrong size.\n\n";

  /*
   *  run ContourPlot::computeIsolinePointsDQ()
   */
  Eigen::MatrixXd isolineDQ = ContourPlot::computeIsolinePointsDQ(G, y0, 5.);
  std::cout << "Output of computeIsolinePointsDQ() is a " << isolineDQ.rows()
            << " by " << isolineDQ.cols() << " matrix";
  if (isolineDQ.cols() > 5)
    std::cout << ", starting with\n" << isolineDQ.leftCols(5) << " ...\n";
  else
    std::cout << ":\n" << isolineDQ << std::endl;

  /*
   *  test ContourPlot::computeIsolinePointsDQ()
   */
  if (isolineDQ.rows() == 2 && isolineDQ.cols() >= 10) {
    // ode45 takes at least 10 steps by default.
    // Check initial state:
    bool passed_computeIsolinePointsDQ = ((y0 - isolineDQ.col(0)).norm() < tol);
    double Gy0 = G(y0);
    for (int i = 1; i < isolineDQ.cols(); i++) {
      // Check that point lies on the G(y0)-isoline:
      double err = std::abs(G(isolineDQ.col(i)) - Gy0);
      // Check that point is not stationary:
      double step = (isolineDQ.col(i) - isolineDQ.col(i - 1)).norm();
      if (step < tol * tol || err > tol) {
        passed_computeIsolinePointsDQ = false;
        break;
      }
    }
    if (passed_computeIsolinePointsDQ)
      std::cout << "Test for computeIsolinePointsDQ() passed!\n\n";
    else
      std::cout << "Test for computeIsolinePointsDQ() failed.\n\n";
  } else
    std::cout << "Test for computeIsolinePointsDQ() failed: wrong size.\n\n";

  return 0;
}
