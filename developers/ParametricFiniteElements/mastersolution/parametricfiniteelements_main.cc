/**
 * @file parametricfiniteelements.cc
 * @brief NPDE homework ParametricFiniteElements code
 * @author Amélie Loher
 * @date 04.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "parametricfiniteelements.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

using namespace ParametricFiniteElements;

int main() {

  // Creating a 2D mesh on the interval [0,1] x [0,Psi(x(0))]
  // nb. of cells
  
  /*unsigned int n = 100;                   
  Eigen::MatrixXd mesh(n+1, 2); 
  
  // Topography function
  auto Psi = [] (double x) {
	  return x * x + 1.0;
  };

  // Nodes are equally spaced in each direction
  for (int i = 0; i <  + 1; i++) {
    mesh[i, 0] = i * (1.0 / n);
	mesh[i, 1] = i * (1.0 / Psi(mesh[i, 0]));
  }
*/
  unsigned int n = 100;
  auto Psi = [](double x) -> double {
	  return x * x + 1.0;
  };

  auto alpha = [](Eigen::Vector2d x) -> double {
	  return 3.0/2.0;
  };

  Eigen::VectorXd mu((n+1)*(n+1));

  mu = geoThermSolve(n, alpha, Psi);

  double val = geoThermSurfInt(n, Psi, mu);

  return 0;
}
