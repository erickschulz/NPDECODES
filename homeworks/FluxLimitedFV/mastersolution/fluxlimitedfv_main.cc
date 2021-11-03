/**
 * @ file
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#define _USE_MATH_DEFINES

#include <Eigen/Core>
#include <cmath>
#include <fstream>
#include <iostream>

#include "fluxlimitedfv.h"

typedef std::numeric_limits<double> dbl;

using namespace FluxLimitedFV;

int main(int /*argc*/, char** /*argv*/) {
  std::cout.precision(dbl::max_digits10);

  /* ADVECTION PROBLEM */
  // Data
  double T = 1.0;
  unsigned int nb_timesteps[8] = {10, 20, 40, 80, 160, 320, 640, 1280};

  // Parameter coefficient
  double beta = 1.0;

  // Flux limiter function
  auto phi = [](double theta) {
    return (abs(theta) + theta) / (1.0 + abs(theta));
  };

  // Initial conditions A and B as lambda functions
  auto mu0_A_fn = [](double x) -> double {
    if (x < 1) {
      return 0.0;
    } else if (x < 2) {
      return pow(sin(0.5 * M_PI * (x - 1)), 2);
    } else {
      return 1.0;
    }
  };
  auto mu0_B_fn = [](double x) -> double {
    if ((1 < x) && (x < 2)) {
      return 1.0;
    } else {
      return 0.0;
    }
  };

  Eigen::VectorXd error_L1_mu0_A = Eigen::VectorXd::Zero(8);
  Eigen::VectorXd error_L1_mu0_B = Eigen::VectorXd::Zero(8);
  double tau, h;
  int N;
  Eigen::VectorXd fluxlimAdvection_sol_A, fluxlimAdvection_sol_B;

  for (int k = 0; k < 8; k++) {
    tau = T / nb_timesteps[k];
    h = 1.2 * tau;
    N = std::round(5.0 / h);

    // Discrete initial conditions
    Eigen::VectorXd mu0_A(N);
    for (int j = 0; j < N; j++) {
      mu0_A(j) = mu0_A_fn(h * j);
    }
    Eigen::VectorXd mu0_B(N);
    for (int j = 0; j < N; j++) {
      mu0_B(j) = mu0_B_fn(h * j);
    }

    fluxlimAdvection_sol_A =
        fluxlimAdvection(beta, mu0_A, h, tau, nb_timesteps[k], phi);
    fluxlimAdvection_sol_B =
        fluxlimAdvection(beta, mu0_B, h, tau, nb_timesteps[k], phi);

    // Computing the L1 errors
    double sum_A = 0.0;
    double sum_B = 0.0;
    for (int j = 0; j < N; j++) {
      sum_A = sum_A +
              std::abs(h * (fluxlimAdvection_sol_A(j) - mu0_A_fn(j * h - T)));
      sum_B = sum_B +
              std::abs(h * (fluxlimAdvection_sol_B(j) - mu0_B_fn(j * h - T)));
    }
    error_L1_mu0_A[k] = sum_A;
    error_L1_mu0_B[k] = sum_B;
  }

  std::cout << "" << std::endl;
  std::cout << "--------------" << std::endl;
  std::cout << "error_L1_mu0_A" << std::endl;
  std::cout << "--------------" << std::endl;
  std::cout << "M";
  std::cout << "\t error" << std::endl;
  for (int k = 0; k < 8; k++) {
    std::cout << nb_timesteps[k];
    std::cout << "\t";
    std::cout << error_L1_mu0_A[k] << std::endl;
  }
  std::cout << "\n" << std::endl;
  std::cout << "--------------" << std::endl;
  std::cout << "error_L1_mu0_B" << std::endl;
  std::cout << "--------------" << std::endl;
  std::cout << "M";
  std::cout << "\t error" << std::endl;
  for (int k = 0; k < 8; k++) {
    std::cout << nb_timesteps[k];
    std::cout << "\t";
    std::cout << error_L1_mu0_B[k] << std::endl;
  }
  std::cout << "" << std::endl;

  // Output results to csv files
  Eigen::VectorXd x_advection(N);
  for (int j = 0; j < N; j++) {
    x_advection(j) = j * h;
  }
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream fluxlimAdvection_sol_A_csv;
  fluxlimAdvection_sol_A_csv.open(CURRENT_BINARY_DIR
                                  "/fluxlimAdvection_sol_A.csv");
  fluxlimAdvection_sol_A_csv << x_advection.transpose().format(CSVFormat)
                             << std::endl;
  fluxlimAdvection_sol_A_csv
      << fluxlimAdvection_sol_A.transpose().format(CSVFormat) << std::endl;
  fluxlimAdvection_sol_A_csv.close();
  std::cout << "Generated " CURRENT_BINARY_DIR "/fluxlimAdvection_sol_A.csv"
            << std::endl;
  std::system("python3 " CURRENT_SOURCE_DIR "/plot_sol.py " CURRENT_BINARY_DIR
              "/fluxlimAdvection_sol_A.csv " CURRENT_BINARY_DIR
              "/fluxlimAdvection_sol_A.eps");

  std::ofstream fluxlimAdvection_sol_B_csv;
  fluxlimAdvection_sol_B_csv.open(CURRENT_BINARY_DIR
                                  "/fluxlimAdvection_sol_B.csv");
  fluxlimAdvection_sol_B_csv << x_advection.transpose().format(CSVFormat)
                             << std::endl;
  fluxlimAdvection_sol_B_csv
      << fluxlimAdvection_sol_B.transpose().format(CSVFormat) << std::endl;
  fluxlimAdvection_sol_B_csv.close();
  std::cout << "Generated " CURRENT_BINARY_DIR "/fluxlimAdvection_sol_B.csv"
            << std::endl;
  std::system("python3 " CURRENT_SOURCE_DIR "/plot_sol.py " CURRENT_BINARY_DIR
              "/fluxlimAdvection_sol_B.csv " CURRENT_BINARY_DIR
              "/fluxlimAdvection_sol_B.eps");

  /* BURGERS FLUX PROBLEM */
  // Discretization parameters
  T = 2.0;
  tau = T / nb_timesteps[7];
  h = 1.2 * tau;
  N = std::round(5.0 / h);

  // Initial conditions A and B as lambda functions
  auto mu0_fn = [](double x) -> double {
    if (x < 2) {
      return -1.0;
    } else if (x < 3) {
      return 1.0;
    } else {
      return -1.0;
    }
  };

  Eigen::VectorXd mu0(N);
  for (int j = 0; j < N; j++) {
    mu0(j) = mu0_fn(h * j);
  }

  Eigen::VectorXd fluxlimBurgers_sol =
      fluxlimBurgers(mu0, h, tau, nb_timesteps[6], phi);

  Eigen::VectorXd x_Burgers(N);
  for (int j = 0; j < N; j++) {
    x_Burgers(j) = j * h;
  }

  std::ofstream fluxlimBurgers_sol_csv;
  fluxlimBurgers_sol_csv.open(CURRENT_BINARY_DIR "/fluxlimBurgers_sol.csv");
  fluxlimBurgers_sol_csv << x_Burgers.transpose().format(CSVFormat)
                         << std::endl;
  fluxlimBurgers_sol_csv << fluxlimBurgers_sol.transpose().format(CSVFormat)
                         << std::endl;
  fluxlimBurgers_sol_csv.close();
  std::cout << "Generated " CURRENT_BINARY_DIR "/fluxlimBurgers_sol.csv"
            << std::endl;
  std::system("python3 " CURRENT_SOURCE_DIR "/plot_sol.py " CURRENT_BINARY_DIR
              "/fluxlimBurgers_sol.csv " CURRENT_BINARY_DIR
              "/fluxlimBurgers_sol.eps");

}  // main
