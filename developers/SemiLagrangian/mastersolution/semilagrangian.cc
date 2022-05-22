/**
 * @ file
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "semilagrangian.h"

namespace SemiLagrangian {

void testfloor_and_division(){
    int M = 80;
    double h = 1.0/80;
    Eigen::Vector2d x (0.504,0.1625);
    std::cout << "j: " << std::floor(x(1)/h) << "(exact: " << x(1)/h << ")" << std::endl;
    std::cout << "j*h: " <<std::floor(x(1)/h)*h  << std::endl;
    std::cout << "x_loc formula from the exercise (direct computation): " << (x(1)-std::floor(x(1)/h)*h)/h << std::endl;
    std::cout << "x_loc fmod: " << std::fmod(x(1),h)/h << std::endl;
    std::cout << "Backward transformation (exercise): "  << std::floor(x(1)/h)*h + ((x(1)-std::floor(x(1)/h)*h)/h)*h << std::endl;
    std::cout << "Backward transformation (fmod): " << std::floor(x(1)/h)*h +  (std::fmod(x(1),h)/h)*h << std::endl;
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

// Auxiliary function of 'semiLagrangeSource'
double evalFEfunction(const Eigen::Vector2d& x, const Eigen::VectorXd& u) {
  int N = u.size();  // assume dofs on boundary already removed
  int root = std::round(std::sqrt(N));
  int M = root + 1;
  double h = 1. / M;

  //compute the location of the square containing x
  int i = std::floor(x(0)/h);
  int j = std::floor(x(1)/h);

  //Check, that x lies in the unit square
  if(i < 0 || i > M-1 || j < 0 || j > M-1){
      std::cerr << "i,j can only be in [0,M-1]" << std::endl;
  }

  //compute local coordinates:
  Eigen::Vector2d x_loc;
  x_loc(0) = (x(0)-i*h)/h;
  x_loc(1) = (x(1)-j*h)/h;

  //Vector of local coefficients:
  Eigen::Vector4d u_loc;

  //Check for boundary dofs and extract correct components of u:
  //Recall clockwise ordering of local dofs, starting in bottom left corner.
  u_loc(0) = (i==0 ||j==0 ) ? 0.0: u((M-1)*(j-1) + (i-1));
  u_loc(1) = (i==(M-1) ||j ==0) ? 0.0: u((M-1)*(j-1) + i);
  u_loc(2) = (i==(M-1) || j ==(M-1)) ? 0.0: u((M-1)*j + i);
  u_loc(3)=  (i==0 || j ==(M-1) )? 0.0: u((M-1)*j  + i-1);


  //evaluate using reference shape functions:
  return u_loc(0) * (1. - x_loc(0)) * (1. - x_loc(1)) +
         u_loc(1) * (1. - x_loc(1) ) * x_loc(0)  +
         u_loc(2) * x_loc(0)  * x_loc(1)  +
         u_loc(3) * (1. - x_loc(0) ) * x_loc(1);

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
  //auto velocity = [](const Eigen::Vector2d& x) {
  //  return (Eigen::Vector2d() << -(x(1) - 0.5), x(0) - 0.5).finished();
  //};
  auto velocity = [](const Eigen::Vector2d& x){
    return Eigen::Vector2d(-x(1),x(0));
  };

  double tau = T / K;  // timestep
  double t = 0;        // initial time
  for (int i = 0; i < K; ++i) {
    Eigen::VectorXd rhs = semiLagrangeSource(u, tau,velocity);
    u = rhs * (M * M);  // inverse of diagonal matrix
  }

  return u;
}
/* SAM_LISTING_END_3 */

} //namespace SemiLagrangian
