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

TEST(SemiLagrangian,evalFEfunction){
  int M=4;
  int N=(M-1)*(M-1); //9
  Eigen::VectorXd u(9);
  u << 1,1,1,1,2,1,1,1,1;

  //check boundary elements:
  EXPECT_NEAR(evalFEfunction(Eigen::Vector2d(0.25*0.7,0.6), u), 0.7, 1E-5);
  EXPECT_NEAR(evalFEfunction(Eigen::Vector2d(0.25*0.7,0.3), u), 0.7, 1E-5);

  EXPECT_NEAR(evalFEfunction(Eigen::Vector2d(0.6,0.25*0.7), u), 0.7, 1E-5);
  EXPECT_NEAR(evalFEfunction(Eigen::Vector2d(0.3,0.25*0.7), u), 0.7, 1E-5);

  EXPECT_NEAR(evalFEfunction(Eigen::Vector2d(1.0-0.25*0.7,0.6), u), 0.7, 1E-5);
  EXPECT_NEAR(evalFEfunction(Eigen::Vector2d(1.0-0.25*0.7,0.3), u), 0.7, 1E-5);
  
  EXPECT_NEAR(evalFEfunction(Eigen::Vector2d(0.6,1.0-0.25*0.7), u), 0.7, 1E-5);
  EXPECT_NEAR(evalFEfunction(Eigen::Vector2d(0.3,1.0-0.25*0.7), u), 0.7, 1E-5);

  //check one of the central elements
  EXPECT_NEAR(evalFEfunction(Eigen::Vector2d(0.4375,0.3125),u), 1.1875,1E-5 );
}

TEST(SemiLagrangian,semiLagrangeSource){
  int M=4;
  int N=(M-1)*(M-1); //9
  double h = 0.25;
  Eigen::VectorXd u(9);
  u << 1,1,1,1,2,1,1,1,1;

  auto v = [](Eigen::Vector2d x){ return Eigen::Vector2d(-1,-2);};
  double tau=0.15;
  Eigen::MatrixXd grid = findGrid(4);

  //reference source
  Eigen::VectorXd f_ref(9);
  f_ref << 1.48,1.32,0.4,0.8,0.8,0.32,0.0,0.0,0.0;
  f_ref = f_ref * 0.0625;
  
  //computed source
  Eigen::VectorXd f_comp = semiLagrangeSource(u,tau,v);

  EXPECT_NEAR((f_ref-f_comp).norm(),0.0,1E-6);
}


}  // namespace SemiLagrangian::test
