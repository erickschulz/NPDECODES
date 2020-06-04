// * Demonstration code for course Numerical Methods for Partial Differential
// * Equations Author: R. Hiptmair, SAM, ETH Zurich
// * Date: May 2020
// Adaptped from MATLAB codes in
// https://svn.id.ethz.ch/sam/Numcourses/rw/matlab/FiniteVolumesNPDE

#ifndef NUMEXP_RUNNER_H
#define NUMEXP_RUNNER_H

#include <cmath>
#include <fstream>
#include <iostream>
/*
 * Functions defining initial data, commonly denoted by u_0
 * BOX: characteristic function of the interval [0,1]
 * WEDGE: tent function with height 1 over supported in [0,1]
 * BUMP: 1-cos^2 bump supported in [0,1]
 * see Exp. 8.8.4.1 from lecture document
 */
static auto box = [](double x) {
  return ((x >= 0.0) && (x <= 1.0)) ? 1.0 : 0.0;
};
static auto wedge = [](double x) {
  return ((x >= 0.0) && (x <= 1.0)) ? (1.0 - 2 * std::abs(x - 0.5)) : 0.0;
};
static auto bump = [](double x) {
  return ((x >= 0.0) && (x <= 1.0)) ? (1.0 - std::pow(std::cos(M_PI * x), 2))
                                    : 0.0;
};

// Flux function for Burgers equations
static auto fb = [](double u) { return 0.5 * u * u; };

// Define numerical fluxes for Burgers equation
// Note: flux function for Burgers ewquation need not be captured
// 1. Local Lax-Friedrich/Rusanov flux for Burgers equation
static auto nfn_lf_burger = [](double v, double w) {
  return 0.5 * (fb(v) + fb(w)) -
         0.5 * std::max(std::abs(v), std::abs(w)) * (w - v);
};
// 2. Godunov flux for Burgers equation
static auto nfn_god_burger = [](double v, double w) {
  if (v < w)
    return fb(std::min(std::max(v, 0.0), w));
  return std::max(fb(v), fb(w));
};

/*
 * @brief driver function for solving 1D scalar conservation laws by means of
 * conservative finite volume methods
 *
 * @tparam EVLFUNCTION functor type
 * @param evl functor object realizing fully discrete evolution
 * @param T final time
 * @param a left boundary of computational spatial domain
 * @param b right boundary of computational spatial domain
 *
 */

template <typename EVLFUNCTION>
void consform_compute(EVLFUNCTION &&evl, std::string filename, double T = 4.0,
                      double a = -1.0, double b = 5.0) {
  std::cout << "Running driver for discrete evolution in conservation form"
            << std::endl;
  // Number of mesh cells in [0,1] for numerical experiments
  std::vector<unsigned int> Nunit_vals{20, 40, 80, 160};

  // Formatted output to file
  const Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols,
                                  ", ", "\n");
  std::ofstream result_file(filename.c_str());
  // Write paramters for numerical experiment
  result_file << "fn = " << filename << ", a = " << a << ", b = " << b
              << ", T = " << T << std::endl;

  // Carry out computations on a sequence of meshes
  for (const unsigned int n : Nunit_vals) {
    // Number of spatial cells inside [a,b]
    const unsigned int N = (unsigned int)((b - a) * n);
    // Progress information
    std::cout << "FV solve on [a,b] = [" << a << ',' << b << "], N = " << N
              << ", T = " << T << std::endl;
    // Invoke actual simulation, which returns a vector of cell averages
    // corresponding to the finite-volume approximation at final time
    Eigen::VectorXd u_final = evl(a, b, N, T);
    // Output another line of data: simulation parameters first, then the
    // simulation result as a row of comma-separated values
    result_file << (u_final.transpose()).format(CSVFormat) << std::endl;
  }
  result_file.close();
  std::cout << "Finished, generated " << filename << std::endl;
  // The data generated by this function can be visualized by means of the
  // Python script consformvis.py
}

#endif
