///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): R. Hiptmair
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <functional>
#include <iomanip>

#include "odeintadapt.h"

int main() {
  using State_t = double;
  using DiscEvolOp = std::function<State_t(double, State_t)>;

  // Differential equation $\dot{y} = y^2$
  auto f = [](double x) { return pow(x, 2); };
  auto norm = [](double x) { return fabs(x); };

  // explicit euler (order 1)
  DiscEvolOp psilow = [&](double h, State_t y) { return y + h * f(y); };

  // explicit trapezoidal (order 2)
  DiscEvolOp psihigh = [&](double h, State_t y) {
    double k1 = f(y);
    double k2 = f(y + h * k1);
    return y + (h / 2.) * (k1 + k2);
  };

  double y0 = 0.5;  // Initial value
  std::vector<std::pair<double, double>> states =
      odeintadapt(psilow, psihigh, y0, 1.9, 0.2, 1e-2, 1e-2, 1e-4, norm);

  // Output of result
  std::cout << std::setw(16) << "t" << std::setw(16) << "y_k" << std::setw(16)
            << "y(t_k)" << std::endl;
  for (auto sol : states) {
    const double tk = sol.first;
    const double soltk = 1.0 / (2.0 - tk);
    std::cout << std::setw(16) << tk << std::setw(16) << sol.second
              << std::setw(16) << soltk << std::endl;
  }
}
