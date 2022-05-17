#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
#include <limits>

using numeric_t = double;
using sparseMatrix_t = Eigen::SparseMatrix<numeric_t>;
using matrix_t = Eigen::Matrix<numeric_t, Eigen::Dynamic, Eigen::Dynamic>;
using vector_t = Eigen::Matrix<numeric_t, Eigen::Dynamic, 1>;
using coord_t = Eigen::Matrix<numeric_t, 2, 1>;

/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR_V, typename FUNCTOR_U0>
double solveTransport(const coord_t& x, int K, double t, FUNCTOR_V&& v,
                      FUNCTOR_U0&& u0) {
  double tau = t / K;  // timestep
  coord_t y = x;       // starting point

  for (int i = 0; i < K; ++i) {
    // A single step of Heun's method, expressed by RK increments
    coord_t k_1 = v(t, x);
    coord_t k_2 = v(t - 2. / 3. * tau, x - 2. / 3. * tau * k_1);
    y -= tau / 4. * k_1 + 3. / 4. * tau * k_2;
  }

  // Test whether trajectory has hit inflow boundary
  if (y(0) >= 0. && y(0) <= 1. && y(1) >= 0. && y(1) <= 1.) {
    return u0(y);
  } else {
    return 0.;
  }
}
/* SAM_LISTING_END_1 */

// Auxiliary function of 'semiLagrangeSource'
matrix_t findGrid(int M) {
  matrix_t grid(2, (M - 1) * (M - 1));

  double h = 1. / M;
  double x1 = h;

  for (int i = 0; i < M - 1; ++i) {
    double x0 = h;
    for (int j = 0; j < M - 1; ++j) {
      coord_t x;
      x << x0, x1;
      grid.col(i * (M - 1) + j) = x;
      x0 += h;
    }
    x1 += h;
  }

  return grid;
}

// Auxiliary function of 'semiLagrangeSource'
double evalFEfunction(const coord_t& x, const vector_t& u) {
  int N = u.size();  // assume dofs on boundary already removed
  int root = std::round(std::sqrt(N));
  int M = root + 1;
  double h = 1. / M;

  // Restore 0-dofs on boundary:
  vector_t u_all = vector_t::Zero((M + 1) * (M + 1));
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
  vector_t coeffs(4);
  for (int i = 0; i < 4; ++i) {
    coeffs(i) = u_all(indices(i));
  }

  // Local coordinates in grid square:
  coord_t x_loc;
  x_loc << std::fmod(x(0), h), std::fmod(x(1), h);

  // Interpolation:
  return coeffs(0) * (1. - x_loc(0) / h) * (1. - x_loc(1) / h) +
         coeffs(1) * (1. - x_loc(1) / h) * x_loc(0) / h +
         coeffs(2) * x_loc(0) / h * x_loc(1) / h +
         coeffs(3) * (1. - x_loc(0) / h) * x_loc(1) / h;
}

/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
vector_t semiLagrangeSource(const vector_t& u_old, double tau, double t,
                            FUNCTOR&& velocity) {
  // Note: components of coefficient vectors are associated
  // with interior nodes only
  int N = u_old.size();  // assume dofs on boundary already removed
  // Extract number of cells in one direction
  int root = std::round(std::sqrt(N));
  if (N != root * root) {
    std::cerr << "The number of dofs should be a perfect square!" << std::endl;
  }
  int M = root + 1;

  matrix_t grid = findGrid(M);
  vector_t f(N);
  for (int i = 0; i < grid.cols(); ++i) {
    // Find grid point corresponding to a degree of freedom
    coord_t x = grid.col(i);
    // Determine location of advected gridpoint
    coord_t y = x - tau * velocity(t, x);
    // Test whether it still lies inside the computational domain
    if (y(0) >= 0. && y(0) <= 1. && y(1) >= 0. && y(1) <= 1.) {
      // Evaluate finite element function from previous timestep
      // at preimage of gridpoint under flow
      f(i) = evalFEfunction(y, u_old);
    } else {
      // Zero, if advected point outside domain
      f(i) = 0.;
    }
  }
  // Finally scale with $h^{-2}$
  return f / (M * M);  // * 1 (from $[0,1]^2$) * 4 (from no. of adjacent
                       // squares) / 4 (from no. of vertices of square)
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
vector_t semiLagrangePureTransport(int M, int K, double T) {
  int N = (M - 1) * (M - 1);  // internal dofs
                              /*
                                  sparseMatrix_t mat(N, N);
                                  for(int i=0; i<N; ++i) {
                                      mat.insert(i,i) = 1. / (M*M); // * 1 (from $[0,1]^2$) * 4 (from no. of
                                 adjacent squares) / 4 (from no. of vertices of square)
                                  }
                              */
  // Coordinates of nodes of the grid
  matrix_t grid = findGrid(M);
  vector_t u(N);
  // Interpolate intial data
  for (int i = 0; i < grid.cols(); ++i) {
    // Position of grid point with index i
    coord_t x = grid.col(i);
    x(0) -= 0.25;
    x(1) -= 0.5;
    if (x.norm() < 0.25) {
      u(i) = std::pow(std::cos(2. * M_PI * x.norm()), 2);
    } else {
      u(i) = 0.;
    }
  }

  // Lambda function giving velocity
  auto velocity = [](double t, const coord_t& x) {
    return (coord_t() << -(x(1) - 0.5), x(0) - 0.5).finished();
  };

  double tau = T / K;  // timestep
  double t = 0;        // initial time
  for (int i = 0; i < K; ++i) {
    t += tau;
    vector_t rhs = semiLagrangeSource(u, tau, t, velocity);
    u = rhs * (M * M);  // inverse of diagonal matrix
  }

  return u;
}
/* SAM_LISTING_END_3 */

int main() {
  double T = M_PI / 2.;

  auto v = [](double t, const coord_t& x) {
    return (coord_t() << -x(1), x(0)).finished();
  };
  auto u0 = [](const coord_t& x) {
    coord_t x0 = x;
    x0(0) -= 0.25;
    x0(1) -= 0.5;
    if (x0.norm() < 0.25) {
      return std::pow(std::cos(2. * M_PI * x0.norm()), 2);
    } else {
      return 0.;
    }
  };

  std::cout << "M"
            << "\t"
            << "K"
            << "\t"
            << "error" << std::endl;

  for (int M = 10; M <= 640; M *= 2) {
    int N = (M - 1) * (M - 1);
    matrix_t grid = findGrid(M);

    for (int K = 10; K <= 640; K *= 2) {
      vector_t u = semiLagrangePureTransport(M, K, T);

      vector_t u_ex(N);
      for (int i = 0; i < grid.cols(); ++i) {
        coord_t x = grid.col(i);
        u_ex(i) = solveTransport(x, K, T, v, u0);
      }

      double err = (u - u_ex).cwiseAbs().maxCoeff();
      std::cout << M << "\t" << K << "\t" << err << std::endl;
    }
  }
}
