/** @file
 * @brief functions for estimating order of convergence
 * @author Ralf Hiptmair
 * @date March 2019
 */

#include "estimated_order_convergence.h"
#include <mgl2/mgl.h>
#include <Eigen/QR>
#include <cassert>

namespace rate_of_convergence {

// Linear fitting of data passed in vectors \texttt{x} and \texttt{y} of equal
// length.
// Returns the 2-vector \Blue{$[\beta_{\ast},\alpha_{\ast}]^{\top}$}, \emph{cf.}
// \eqref{eq:linfit}.
Eigen::Vector2d linearFit(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
  assert(x.rows() == y.rows());
  Eigen::Matrix<double, Eigen::Dynamic, 2> X(x.rows(), 2);
  // Set up matrix of overdetermined system of equations
  X.col(0) = Eigen::VectorXd::Constant(x.rows(), 1);
  X.col(1) = x;
  // Solve least squares problem by QR-decomposition
  return X.fullPivHouseholderQr().solve(y);
}

double eoc(const Eigen::VectorXd& N, const Eigen::VectorXd& err,
           unsigned fromindex, std::string filename) {
  // The argument \texttt{N} has to pass a \emph{sorted} vector of length
  // \Blue{$L>1$} of problem size parameter values, whereas the
  // \Blue{$L$}-vector \texttt{err} contains the corresponding error norms. The
  // argument \texttt{fromindex}\Blue{$\in \{1,\ldots,L-1\}$} restricts the
  // relevant data to \Blue{$\mathtt{fromindex},\ldots,L\}$} in order to
  // suppress the impact of possible pre-asymptotic behavior
  // Returns the estimated rate of convergence.
  const unsigned dim = N.size();
  // Consistency check for arguments
  assert(dim > 1);
  assert(fromindex + 1 < dim);
  assert(err.size() == dim);

  // check if data is proper, that is, positive
  assert(
      std::none_of(N.data(), N.data() + dim, [](double d) { return d <= 0; }));
  assert(std::none_of(err.data(), err.data() + dim,
                      [](double d) { return d <= 0; }));

  // check no two elements in N are equal, sorting assumed!
  assert(std::is_sorted(N.data(), N.data() + dim));
  assert(N.data() + dim == std::adjacent_find(N.data(), N.data() + dim));

  // truncate preasymptotic behavior if desired:
  const unsigned newdim = dim - fromindex;

  // compute $\log(N)$ and $\log(\mathtt{err})$ componentwise
  auto logfun = [](double d) { return std::log(d); };
  Eigen::VectorXd Nlog(newdim), errlog(newdim);
  std::transform(N.data() + fromindex, N.data() + dim, Nlog.data(), logfun);
  std::transform(err.data() + fromindex, err.data() + dim, errlog.data(),
                 logfun);

  // perform \com{linear regression}, aka least squares fitting to a line.
  // \texttt{linearFit} returns the coefficients of the linear polynomial, the
  // second of which is its slope
  Eigen::Vector2d polyfit = linearFit(Nlog, errlog);
  double alpha = -polyfit[1];
  double offset = std::exp(polyfit[0]);

  // evaluate the polynomial
  Eigen::VectorXd polyval(dim);
  std::transform(N.data(), N.data() + dim, polyval.data(),
                 [&](double d) { return offset * std::pow(d, -alpha); });
  std::string lbl = "O(N^{-" + std::to_string(alpha) + "})";
  plotConvergence(N, err, polyval, lbl, filename);

  return alpha;
}

void plotConvergence(const Eigen::VectorXd& N, const Eigen::VectorXd& err,
                     const Eigen::VectorXd& fit, const std::string& label,
                     const std::string& filename) {
  // Data in mgl compatible format
  mglData N_d(N.data(), N.size()), err_d(err.data(), err.size()),
      fit_d(fit.data(), fit.size());
  mglGraph gr;          // Graph object
  gr.SetFontSizePT(8);  // Smaller font size than default
  gr.SetRanges(N_d.Minimal(), N_d.Maximal(), err_d.Minimal(), err_d.Maximal());
  gr.SetFunc("lg(x)", "lg(y)");     // Set loglog
  gr.AddTick('y', 0.01, "0.01");    // Somehow the loglog breaks the automatic
  gr.AddTick('y', 0.1, "0.1");      // axis labels (for small values), so we set
  gr.AddTick('y', 0.5, "0.5");      // them manually here
  gr.Label('x', "Problem size N");  // x label
  gr.Label('y', "Error norm");      // y label
  gr.Grid("", "h");                 // Activate thin grey grid
  gr.Axis();                        // Activate axis display
  gr.Plot(N_d, err_d, " r*");       // Plot the data (whitespace = no line)
  gr.AddLegend("Error", " r*");     // Label for the data
  gr.Plot(N_d, fit_d, "b0");        // Plot the fit as line
  gr.AddLegend(label.c_str(), "b0");  // Label for the fit
  gr.Legend();                        // Add legend
  gr.WriteEPS(filename.c_str());      // Save figure
}

}  // namespace rate_of_convergence
