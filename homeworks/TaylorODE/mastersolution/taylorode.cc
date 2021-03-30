/**
 * @file taylorode.cc
 * @brief NPDE homework TaylorODE
 * @author ?, Philippe Peter
 * @date 24.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "taylorode.h"

#include <Eigen/Core>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace TaylorODE {

Eigen::Vector2d PredPreyModel::f(const Eigen::Vector2d& y) const {
  return {(alpha1_ - beta1_ * y(1)) * y(0), (beta2_ * y(0) - alpha2_) * y(1)};
}

Eigen::Vector2d PredPreyModel::df(const Eigen::Vector2d& y,
                                  const Eigen::Vector2d& z) const {
  Eigen::Matrix2d Df;
  Df << alpha1_ - beta1_ * y(1), -beta1_ * y(0), beta2_ * y(1),
      -alpha2_ + beta2_ * y(0);
  return Df * z;
}

Eigen::Vector2d PredPreyModel::d2f(const Eigen::Vector2d& y,
                                   const Eigen::Vector2d& z) const {
  Eigen::Matrix2d H1, H2;
  H1 << 0, -beta1_, -beta1_, 0;
  H2 << 0, beta2_, beta2_, 0;
  return {z.transpose() * H1 * z, z.transpose() * H2 * z};
}

std::vector<Eigen::Vector2d> SolvePredPreyTaylor(const PredPreyModel& model,
                                                 double T,
                                                 const Eigen::Vector2d& y0,
                                                 unsigned int M) {
  std::vector<Eigen::Vector2d> res;
  res.reserve(M + 1);

  Eigen::Vector2d y = y0;
  double h = T / M;
  res.push_back(y);

  for (unsigned int k = 0; k < M; ++k) {
    // evaluate terms for taylor step.
    auto fy = model.f(y);              // f(y)
    auto dfyfy = model.df(y, fy);      // df(y) * f(y)
    auto df2yfy = model.df(y, dfyfy);  // df(y)^2 * f(y)
    auto d2fyfy = model.d2f(y, fy);    // d2f(y)(fy,fy)

    // evaluate taylor expansion to compute update
    y = y + h * fy + 0.5 * h * h * dfyfy +
        1.0 / 6.0 * h * h * h * (df2yfy + d2fyfy);

    // save new state:
    res.push_back(y);
  }
  return res;
}

double TestCvgTaylorMethod() {
  // initialize parameters for the model:
  double T = 10;               // final time
  Eigen::Vector2d y0(100, 5);  // initial condition
  Eigen::Vector2d yex(0.319465882659820,
                      9.730809352326228);  // reference solution
  double alpha1 = 3.0;
  double alpha2 = 2.0;
  double beta1 = 0.1;
  double beta2 = 0.1;
  PredPreyModel model(alpha1, alpha2, beta1, beta2);

  // Initialize parameters for the convergence study
  unsigned int M0 = 128;    // Minimum number of timesteps
  unsigned int numRef = 8;  // Number of refinements

  // Convergence study
  Eigen::ArrayXd error(numRef);
  Eigen::ArrayXd M(numRef);
  for (unsigned int i = 0; i < numRef; ++i) {
    M(i) = std::pow(2, i) * M0;
    auto res = SolvePredPreyTaylor(model, T, y0, M(i));
    error(i) = (res.back() - yex).norm();
  }

  PrintErrorTable(M, error);

  // calculate linear regression line: log(error) ~ c0 + c1*log(M)
  Eigen::MatrixXd A(numRef, 2);
  A.col(0) = Eigen::VectorXd::Ones(numRef);
  A.col(1) = M.log();
  Eigen::VectorXd logError = error.log();
  Eigen::Vector2d coeffs = A.householderQr().solve(logError);

  // estimated convergence rate: -c1
  return -coeffs(1);
}

void PrintErrorTable(const Eigen::ArrayXd& M, const Eigen::ArrayXd& error) {
  std::cout << std::setw(15) << "M" << std::setw(15) << "error" << std::setw(15)
            << "rate" << std::endl;

  for (unsigned int i = 0; i < M.size(); ++i) {
    std::cout << std::setw(15) << M(i) << std::setw(15) << error(i);
    if (i > 0) {
      std::cout << std::setw(15) << std::log2(error(i - 1) / error(i));
    }
    std::cout << std::endl;
  }
}

}  // namespace TaylorODE