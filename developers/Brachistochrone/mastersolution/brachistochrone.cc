/**
 * @file XXX.cc
 * @brief NPDE homework XXX code
 * @author 
 * @date 
 * @copyright Developed at SAM, ETH Zurich
 */

#include "brachistochrone.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>

namespace Brachistochrone {
    Eigen::VectorXd coeff_sigma(
        const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots
    )
    {
        int numCells = knots.cols()-1;
        Eigen::VectorXd sigma(numCells);
        double h = 1./(numCells);
        for(int i=0; i<numCells; ++i){
            sigma(i) = 1./(std::sqrt(-.5*(knots(1,i)+knots(1,i+1)))*1./h*(knots.col(i+1)-knots.col(i)).norm());
        }
        //std::cout << "sigma: " << sigma << std::endl;
        return sigma;
    }

    Eigen::VectorXd sourcefn2(
        const Eigen::Matrix<double,2,Eigen::Dynamic> &knots
    )
    {
        //std::cout << knots << std::endl;
        int numCells = knots.cols()-1;
        Eigen::VectorXd f(numCells);
        double h = 1./(numCells);
        for(int i=0; i<numCells; ++i){
            Eigen::VectorXd uhp = 1./h*(knots.col(i+1)-knots.col(i));
            double uh2 = .5*(knots(1,i)+knots(1,i+1));
            //std::cout << uhp << std::endl;
            //f(i) = 1./h*(knots.col(i+1)-knots.col(i)).norm()/(2.*std::sqrt(-.5*(knots(1,i)+knots(1,i+1)))*.5*(knots(1,i)+knots(1,i+1)));
            f(i) = uhp.norm()/(2.*std::sqrt(-uh2)*uh2);
        }
        //std::cout << "f:" << f << std::endl;
        return f;
    }

    Eigen::SparseMatrix<double> matR(
        const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots
    )
    {
        int M = knots.cols();

        double h = 1./(M-1);
        Eigen::SparseMatrix<double> R(M-2,M-2);
        R.reserve(Eigen::RowVectorXi::Constant(M-2,3));
        Eigen::VectorXd sigma = coeff_sigma(knots);
        R.insert(0,0) = (sigma(0)+sigma(1))/h;
        R.insert(0,1) = -sigma(1)/h;
        for(int i=1; i<M-3; ++i){
            R.insert(i,i-1) = -sigma(i)/h;
            R.insert(i,i)   =  (sigma(i)+sigma(i+1))/h;
            R.insert(i,i+1) = -sigma(i+1)/h;
        }
        R.insert(M-3,M-4) = -sigma(M-3)/h;
        R.insert(M-3,M-3) = (sigma(M-3)+sigma(M-2))/h;
        

        //R.insert(0,0) = 2.;
        //R.insert(0,1) = -1;
        //for(int i=1; i<M-3; ++i){
        //    R.insert(i,i-1) = -1;
        //    R.insert(i,i)   =  2.;
        //    R.insert(i,i+1) = -1.;
        //}
        //R.insert(M-3,M-4) = -1.;
        //R.insert(M-3,M-3) = 2.;
        R.makeCompressed();

        return R;
    }

    Eigen::VectorXd compute_rhs(
        const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots,
        Eigen::Vector2d b
    )
    {
        int numCells = knots.cols()-1;
        double h = 1./numCells;
        Eigen::VectorXd sigma = coeff_sigma(knots);
        Eigen::VectorXd f = sourcefn2(knots);

        Eigen::VectorXd phi(2*(numCells-1));
        phi.setZero();
        phi(numCells-2) = sigma(numCells-1)*b(0)/h;

        phi(numCells-1) = sigma(0)*(-0.)/h+.5*h*(f(0)+f(1));
        for(int i=1;i<numCells-2;++i){
            phi(numCells-1+i) = .5*h*(f(i)+f(i+1));
        }
        phi(2*numCells-3) = sigma(numCells-1)*b(1)/h + .5*h*(f(numCells-2)+f(numCells-1));

        return phi;

    }
    
}  // namespace Brachistochrone
