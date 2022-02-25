/**
 * @ file transformedconslaw.cc
 * @ brief NPDE homework about conservation law with non-linear density
 * @ author Ralf Hiptmair, Oliver Rietmann
 * @ date July 2021
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "transformedconslaw.h"

#include <cmath>

namespace TRFCL {

// Use functions rho and g so that f(u)=0.5*u^2.

// Function rho
double NonStdCauchyProblemCL::rho(double z) const { return std::exp(z); }
// Derivative of function rho
double NonStdCauchyProblemCL::drho(double z) const { return rho(z); }
// Function g
double NonStdCauchyProblemCL::g(double z) const {
  return 0.5 * std::exp(2.0 * z);
}  // return 0.5 * z * z;
// Derivative of function g
double NonStdCauchyProblemCL::dg(double z) const {
  return std::exp(2.0 * z);
}  // return z;
// Finite interval containing "interesting" parts of solution
std::pair<double, double> NonStdCauchyProblemCL::domain() const {
  return {-3.0, 3.0};
}
// Final time
double NonStdCauchyProblemCL::T() const { return 1.0; }
// Initial data
double NonStdCauchyProblemCL::z0(double x) const {
  return (std::abs(x) > 1.0) ? 0.0 : std::cos(M_PI * x / 2.0);
}

}  // namespace TRFCL
