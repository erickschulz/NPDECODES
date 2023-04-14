/**
 * @file XXX_test.cc
 * @brief NPDE homework XXX code
 * @author 
 * @date 
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../brachistochrone.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Brachistochrone::test {

    TEST(Brachistochrone, coeff_sigma) {
        Eigen::Matrix<double,2,3> knots = (Eigen::Matrix<double,2,3>() << 0,-1.,-2.,0.,-.5,-1.).finished();
        Eigen::VectorXd sigma = Brachistochrone::coeff_sigma(knots);

        double tol = 1e-2;
        EXPECT_NEAR(sigma(0), 1./std::sqrt(1.25), tol);
        EXPECT_NEAR(sigma(1), 1./std::sqrt(3.75), tol);
    }

    TEST(Brachistochrone, sourcefn2) {
        Eigen::Matrix<double,2,3> knots = (Eigen::Matrix<double,2,3>() << 0,-1.,-2.,0.,-.5,-1.).finished();
        Eigen::VectorXd f = Brachistochrone::sourcefn2(knots);

        double tol = 1e-2;
        EXPECT_NEAR(f(0), std::sqrt(5.)/(-2.*std::sqrt(.25)*.25), tol);
        EXPECT_NEAR(f(1), std::sqrt(5.)/(-2.*std::sqrt(.75)*.75), tol);
    }

    TEST(Brachistochrone, matR) {
        Eigen::Matrix<double,2,6> knots = (Eigen::Matrix<double,2,6>() << 0,-1.,-2.,-3.,-4.,-5.,0.,-.5,-1.,-1.5,-2.,-2.5).finished();
        Eigen::SparseMatrix<double> R = Brachistochrone::matR(knots);

        EXPECT_EQ(R.nonZeros(),10);
        Eigen::MatrixXd Rd=R.toDense();
        Eigen::VectorXd sigma = Brachistochrone::coeff_sigma(knots);
        int M = knots.cols();
        double h = 1./(M-1);
        double tol = 1e-8;
        for(int i=0;i<M-2;++i){
            EXPECT_NEAR(Rd(i,i),(sigma(i)+sigma(i+1))/h,tol);
            if(i>0) EXPECT_NEAR(Rd(i,i-1),-sigma(i)/h,tol);
            if(i<M-3) EXPECT_NEAR(Rd(i,i+1),-sigma(i+1)/h,tol);
        }
    }

    TEST(Brachistochrone, compute_rhs) {
        Eigen::Matrix<double,2,6> knots = (Eigen::Matrix<double,2,6>() << 0,-1.,-2.,-3.,-4.,-5.,0.,-.5,-1.,-1.5,-2.,-2.5).finished();
        Eigen::SparseMatrix<double> R = Brachistochrone::matR(knots);

        EXPECT_EQ(R.nonZeros(),10);
        Eigen::MatrixXd Rd=R.toDense();
        Eigen::Vector2d b = {1,-2.};
        Eigen::VectorXd sigma = Brachistochrone::coeff_sigma(knots);
        Eigen::VectorXd f = Brachistochrone::sourcefn2(knots);
        Eigen::VectorXd rhs = Brachistochrone::compute_rhs(knots,b);
        int M = knots.cols();
        double tol = 1e-8;
        for(int i=0;i<M-3;++i)
            EXPECT_NEAR(rhs(i),0.,tol);
        EXPECT_NEAR(rhs(M-3), -b(0)*sigma(M-2)*(M-1),tol);
        for(int i=0;i<M-3;++i)
            EXPECT_NEAR(rhs(M-2+i),.5/(M-1)*(f(i)+f(i+1)),tol);
        EXPECT_NEAR(rhs(2*M-5), -sigma(M-2)*b(1)*(M-1)+.5/(M-1)*(f(M-3)+f(M-2)),tol);

    }

    TEST(Brachistochrone, brachistochrone) {
        Eigen::Vector2d b = {M_PI,-2};


        int M = 50;
        int atol = 1e-8;
        int rtol = 1e-8;
        unsigned itmax = 100;

        Eigen::MatrixXd knots(2,M+1);
        for(int i=0; i<M+1;++i) knots(0,i) = (M_PI*i)/M-std::sin((M_PI*i)/M);
        for(int i=0; i<M+1;++i) knots(1,i) = std::cos((M_PI*i)/M)-1.;



        Brachistochrone::brachistochrone(M,b,atol,rtol,itmax);

        //std::cout << knots << std::endl;


    }

}  // namespace Brachistochrone::test
