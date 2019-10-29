/**
 * @file discontinuousgalerkin1d.h
 * @brief NPDE homework "DiscontinuousGalerkin1D" code
 * @author Oliver Rietmann
 * @date 22.05.2019
 * @copyright Developed at ETH Zurich
 */

#include <cmath>
#include <utility>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace DiscontinuousGalerkin1D {

/**
 * @brief Returns the mass matrix B according to the exercise.
 * @param Ml number of negative spacial nodes
 * @param Mr number of positive spacial nodes
 * @param h equidistant spacial mesh-width
 * @return diagonal square matrix of dimension 2 * (Ml + Mr + 1)
 */
Eigen::SparseMatrix<double> compBmat(int Ml, int Mr, double h);

/**
 * @brief Computes G according to the problem description.
 * @param mu expansion coefficients vector of size 2 * (Ml + Mr + 1)
 * @param f flux function matching std::function<double(double)>
 * @param F numerical flux matching std::function<double(double, double)>
 * @param Ml number of negative spacial nodes
 * @param Mr number of positive spacial nodes
 * @param h equidistant spacial mesh-width
 * @return vector of length 2 * (Ml + Mr + 1) approximating G(mu)
 */
template <typename FUNCTOR, typename NUMFLUX>
Eigen::VectorXd G(const Eigen::VectorXd &mu, FUNCTOR &&f, NUMFLUX &&F, int Ml,
                  int Mr, double h) {
  const int N_half = (Ml + Mr + 1);
  const int N = 2 * N_half;
  Eigen::VectorXd Gvec(N);
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return Gvec;
}

/**
 * @brief Time evolution by an explicit 2nd order Runge-Kutta method.
 * @param mu0 expansion coefficients vector of size 2 * (Ml + Mr + 1)
 * representing the initial data
 * @param f flux function matching std::function<double(double)>
 * @param F numerical flux matching std::function<double(double, double)>
 * @param Ml number of negative spacial nodes
 * @param Mr number of positive spacial nodes
 * @param h equidistant spacial mesh-width
 * @return expansion coefficients vector of length 2 * (Ml + Mr + 1)
 * representing the solution at time T
 */
template <typename FUNCTOR, typename NUMFLUX>
Eigen::VectorXd dgcl(Eigen::VectorXd mu0, FUNCTOR &&f, NUMFLUX &&F, double T,
                     int Ml, int Mr, double h, unsigned int m) {
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return mu0;
}

/**
 * @brief Engquist-Osher numerical flux for the flux function f(u) = u * (1 -
 * u).
 * @param mu0 expansion coefficients vector of size 2 * (Ml + Mr + 1)
 * representing the initial data
 * @param v left state
 * @param w right state
 * @return F(v, w)
 */
double Feo(double v, double w);

/**
 * @brief Time evolution according to dgcl(...) on the spacial interval [-2, 2]
 * with mesh-width h = 0.05, on th time interval [0, 1] with timestep size h
 * / 3. The solution at endtime is written to solution.csv.
 */
void solveTrafficFlow();

}  // namespace DiscontinuousGalerkin1D
