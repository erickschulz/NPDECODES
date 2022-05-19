/**
 * @ file
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "semilagrangian.h"

namespace SemiLagrangian {


// Auxiliary function of 'semiLagrangeSource'
double evalFEfunction(const Eigen::Vector2d& x, const Eigen::VectorXd& u) {
  int N = u.size();  // assume dofs on boundary already removed
  int root = std::round(std::sqrt(N));
  int M = root + 1;
  double h = 1. / M;

  // Restore 0-dofs on boundary:
  Eigen::VectorXd u_all = Eigen::VectorXd::Zero((M + 1) * (M + 1));
  unsigned cnt = 0;
  for (int i = 0; i < u_all.size(); ++i) {
    if ((i > M) && ((i + 1) % (M + 1) != 0) && (i < M * (M + 1)) &&
        (i % (M + 1) != 0)) {
      u_all(i) = u(cnt);
      ++cnt;
    }
  }

  Eigen::Vector2i x_idx;
  x_idx << std::floor(x(0) / h), std::floor(x(1) / h);
  // Index of bottom-left corner of grid square:
  int ii = x_idx(1) * (M + 1) + x_idx(0);
  if (ii < 0 || ii >= (M + 1) * (M + 1)) {
    std::cerr << "ii can only be in [0,M*M+2*M]" << std::endl;
  }

  // Indices of corners of grid square:
  Eigen::Vector4i indices;
  indices << ii, ii + 1, ii + M + 2, ii + M + 1;

  // Coefficients of corners of grid square:
  Eigen::VectorXd coeffs(4);
  for (int i = 0; i < 4; ++i) {
    coeffs(i) = u_all(indices(i));
  }

  // Local coordinates in grid square:
  Eigen::Vector2d x_loc;
  x_loc << std::fmod(x(0), h), std::fmod(x(1), h);

  // Interpolation:
  return coeffs(0) * (1. - x_loc(0) / h) * (1. - x_loc(1) / h) +
         coeffs(1) * (1. - x_loc(1) / h) * x_loc(0) / h +
         coeffs(2) * x_loc(0) / h * x_loc(1) / h +
         coeffs(3) * (1. - x_loc(0) / h) * x_loc(1) / h;
}


// Auxiliary function of 'semiLagrangeSource'
Eigen::MatrixXd findGrid(int M) {
  Eigen::MatrixXd grid(2, (M - 1) * (M - 1));

  double h = 1. / M;
  double x1 = h;

  for (int i = 0; i < M - 1; ++i) {
    double x0 = h;
    for (int j = 0; j < M - 1; ++j) {
      Eigen::Vector2d x;
      x << x0, x1;
      grid.col(i * (M - 1) + j) = x;
      x0 += h;
    }
    x1 += h;
  }

  return grid;
}


/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd semiLagrangePureTransport(int M, int K, double T) {
  int N = (M - 1) * (M - 1);  // internal dofs
                              /*
                                  sparseMatrix_t mat(N, N);
                                  for(int i=0; i<N; ++i) {
                                      mat.insert(i,i) = 1. / (M*M); // * 1 (from $[0,1]^2$) * 4 (from no. of
                                 adjacent squares) / 4 (from no. of vertices of square)
                                  }
                              */
  // Coordinates of nodes of the grid
  Eigen::MatrixXd grid = findGrid(M);
  Eigen::VectorXd u(N);
  // Interpolate intial data
  for (int i = 0; i < grid.cols(); ++i) {
    // Position of grid point with index i
    Eigen::Vector2d x = grid.col(i);
    x(0) -= 0.25;
    x(1) -= 0.5;
    if (x.norm() < 0.25) {
      u(i) = std::pow(std::cos(2. * M_PI * x.norm()), 2);
    } else {
      u(i) = 0.;
    }
  }

  // Lambda function giving velocity
  auto velocity = [](double t, const Eigen::Vector2d& x) {
    return (Eigen::Vector2d() << -(x(1) - 0.5), x(0) - 0.5).finished();
  };

  double tau = T / K;  // timestep
  double t = 0;        // initial time
  for (int i = 0; i < K; ++i) {
    t += tau;
    Eigen::VectorXd rhs = semiLagrangeSource(u, tau, t, velocity);
    u = rhs * (M * M);  // inverse of diagonal matrix
  }

  return u;
}
/* SAM_LISTING_END_3 */

} //namespace SemiLagrangian
