#ifndef CONSFORMEVL_HPP
#define CONSFORMEVL_HPP

#include "ode.h"

#include <Eigen/Dense>

namespace ConsFV {
using namespace Eigen;

// \ref{cpp:fluxdiff} in lecture notes

/* SAM_LISTING_BEGIN_1 */
// arguments: (Finite) state vector \Blue{$\mu$} of cell averages, see
// \eqref{eq:cellavg}
//  Functor \Blue{$F : \mathbb R \times \mathbb R \mapsto \mathbb R$}, 2-point
//  numerical flux
// return value:  Vector with differences of numerical fluxes, which provides
// the right hand side of \eqref{eq:2pcf}
template <typename FunctionF>
VectorXd fluxdiff(const VectorXd& mu, FunctionF F) {
  unsigned n = mu.size();           // length of state vector
  VectorXd fd = VectorXd::Zero(n);  // return vector

  // constant continuation of data for \Blue{$x\leq a$}!
  fd[0] = F(mu[0], mu[1]) - F(mu[0], mu[0]);
  for (unsigned j = 1; j < n - 1; ++j) {
    fd[j] = F(mu[j], mu[j + 1]) - F(mu[j - 1], mu[j]);  // see \eqref{eq:2pcf}
  }
  // constant continuation of data for \Blue{$x\geq b$}!
  fd[n - 1] = F(mu[n - 1], mu[n - 1]) - F(mu[n - 2], mu[n - 1]);
  // Efficient thanks to return value optimization (RVO)
  return fd;
}
/* SAM_LISTING_END_1 */

// \ref{cpp:consformevl} in lecture notes
//
/* SAM_LISTING_BEGIN_2 */
// arguments:
//   Real numbers \Blue{$a, b$}, the boundaries of the interval,
//   unsigned int \Blue{$N$}, the number of cells,
//   Functor \Blue{$u0 : \mathbb R \mapsto \mathbb R$}, initial value,
//   Final time \Blue{$T>0$},
//   Functor \Blue{$F=F(v,w)$} for 2-point numerical flux function.
//
// return value:
//   Vector with cell values at final time \Blue{$T$}
//
// Finite volume discrete evolution in conservation form with 2-point
// flux, see \eqref{eq:2pcf}; Cauchy problem over time \Blue{$[0,T]$}
template <typename FunctionU0, typename FunctionF>
VectorXd consformevl(double a, double b, unsigned N, FunctionU0 u0, double T,
                     FunctionF F) {
  double h = (b - a) / N;  // meshwidth
  // centers of dual cells
  VectorXd x = VectorXd::LinSpaced(N, a + 0.5 * h, b - 0.5 * h);

  // vector \Blue{$\vec{\mubf}_0$} of initial cell averages
  // obtained by point sampling of \Blue{$u_0$} in grid points
  VectorXd mu0 = x.unaryExpr(u0);

  // right hand side function for ode solver
  auto odefun = [&](const VectorXd& mu, VectorXd& dmdt, double t) {
    dmdt = -1. / h * fluxdiff<FunctionF>(mu, F);
  };

  // Method of lines approach, \emph{c.f.} Sect.~\ref{sec:pmol}: timestepping by
  // Boost integrator (adaptive explicit embedded Runge-Kutta method
  // of order 5, see also Def.~\ref{def:rk2})
  double abstol = 1E-8, reltol = 1E-6;  // integration control parameters
  VectorXd t;                           // Returns temporal grid
  MatrixXd MU;                          // Returns state matrix
  std::tie(t, MU) = ode45(odefun, 0, T, mu0, abstol,
                          reltol);  // \Label[line]{cle:ode45}
  // Final state vector is the rightmost column of MU.
  return MU.col(t.size() - 1);
}
/* SAM_LISTING_END_2 */

}  // namespace ConsFV
#endif
