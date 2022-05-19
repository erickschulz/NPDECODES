/**
 * @file semilagrangian_test.cc
 * @brief NPDE homework TEMPLATE MAIN FILE
 * @author
 * @date May 2022
 * @copyright Developed at SAM, ETH Zurich
 */

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <iostream>

#include "../semilagrangian.h"

namespace SemiLagrangian::test {

TEST(SemiLagrangian, solveTransport_constant_velocity){
  double t = 0.5;
  int K = 100; 
  auto u0 = [](const Eigen::Vector2d& x) {
    Eigen::Vector2d x0 = x;
    x0(0) -= 0.25;
    x0(1) -= 0.5;
    if (x0.norm() < 0.25) {
      return std::pow(std::cos(2. * M_PI * x0.norm()), 2);
    } else {
      return 0.;
    }
  };

  auto v = [](const Eigen::Vector2d& x){
    return Eigen::Vector2d(0.5,0.5);
  };
  auto u_exact = [&u0,&v](double t, const Eigen::Vector2d& x){
    Eigen::Vector2d y = x - v(x)*t;
    if(y(0) >= 0 && y(0) <= 1 && y(1) >= 0 && y(1) <= 1){
      return u0(y);
    }
    else{
      return 0.0;
    }
  };

  int M = 100;
  double h = 1.0/M;
  for(int i = 0; i < M; ++i){
    for(int j=0; j < M; ++j){
      Eigen::Vector2d x (i*h,j*h);
      EXPECT_NEAR(u_exact(t,x),solveTransport(x,K,t,v,u0), 1E-5);
    }
  }
}

TEST(SemiLagrangian, solveTransport_rotation){
  double t = 0.5;
  int K = 100; 
  auto u0 = [](const Eigen::Vector2d& x) {
    Eigen::Vector2d x0 = x;
    x0(0) -= 0.25;
    x0(1) -= 0.5;
    if (x0.norm() < 0.25) {
      return std::pow(std::cos(2. * M_PI * x0.norm()), 2);
    } else {
      return 0.;
    }
  };

  auto v = [](const Eigen::Vector2d& x){
    return Eigen::Vector2d(-x(1),x(0));
  };
  auto u_exact = [&u0,&v](double t, const Eigen::Vector2d& x){
    Eigen::Matrix2d R_inv;
    R_inv << std::cos(-t), -std::sin(-t),
             std::sin(-t), std::cos(-t);
    Eigen::Vector2d y = R_inv*x;
    if(y(0) >= 0 && y(0) <= 1 && y(1) >= 0 && y(1) <= 1){
      return u0(y);
    }
    else{
      return 0.0;
    }
  };

  int M = 100;
  double h = 1.0/M;
  for(int i = 0; i < M; ++i){
    for(int j=0; j < M; ++j){
      Eigen::Vector2d x (i*h,j*h);
      EXPECT_NEAR(u_exact(t,x),solveTransport(x,K,t,v,u0), 1E-4);
    }
  }
}

}  // namespace SemiLagrangian::test
