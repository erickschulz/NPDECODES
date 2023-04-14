/**
 * @file XXX.h
 * @brief NPDE homework XXX code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>
#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

namespace Brachistochrone {

    Eigen::VectorXd coeff_sigma(
        const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots
    );

    Eigen::VectorXd sourcefn2(
        const Eigen::Matrix<double,2,Eigen::Dynamic> &knots
    );

    Eigen::SparseMatrix<double> matR(
        const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots
    );

    Eigen::VectorXd compute_rhs(
        const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots,
        Eigen::Vector2d b
    );




    template <typename RECORDER = std::function<void(const Eigen::Matrix<double, 2, Eigen::Dynamic> &)>>
    Eigen::Matrix<double, 2, Eigen::Dynamic> brachistochrone(
        unsigned int M,
        Eigen::Vector2d b,
        double atol, 
        double rtol,
        unsigned int itmax,
        RECORDER &&rec = [](const Eigen::Matrix<double, 2, Eigen::Dynamic> &) -> void{return;}
    )
    {
        Eigen::MatrixXd knots(2,M+1);
        //for(int i=0; i<M+1; ++i) knots(0,i) = (b(0)*i)/M;
        //for(int i=0; i<M+1; ++i) knots(1,i) = (b(1)*i)/M;
        std::cout << b << std::endl;
        for(int i=0; i<M+1;++i) knots(0,i) = (M_PI*i)/M-std::sin((M_PI*i)/M);
        for(int i=0; i<M+1;++i) knots(1,i) = std::cos((M_PI*i)/M)-1.;
        std::cout << knots << std::endl;
        Eigen::MatrixXd knots_prev = knots;


        auto lencomp = [&knots]() { 
            double out = 0.;
            for(int i=0; i<knots.cols()-1; ++i){
              out+=(knots.col(i)-knots.col(i+1)).norm()/std::sqrt(-.5*(knots(1,i)+knots(1,i+1)));
            }
            return out;
         };

        int q = 0;
        double it_err = 10.;
        do{
            Eigen::SparseMatrix<double> R = matR(knots);
            Eigen::VectorXd rhs = compute_rhs(knots,b);
            //std::cout << "R=" << R.toDense() << std::endl;
            //std::cout << "phi=" << rhs << std::endl;
            //std::cout << std::endl;
            //std::cout << rhs << std::endl;
            //Eigen::VectorXd phi(2);
            //phi(0) = rhs(2);
            //phi(1) = rhs(3);

            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(R);

            LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
            Eigen::VectorXd uh = solver.solve(rhs.segment(0,M-1));
            LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
            //for(int i=0; i<M-1; ++i) knots(0,i+1) = uh(i);
            //std::cout << uh <<std::endl;

            solver.compute(R);
            LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
            Eigen::VectorXd uh2 = solver.solve(rhs.segment(M-1,M-1));
            LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
            for(int i=0; i<M-1; ++i) knots(1,i+1) = uh2(i);
            //std::cout << std::endl;

            //std::cout << knots << std::endl;
            std::cout << std::endl;

            it_err = (knots-knots_prev).colwise().norm().sum()/knots.colwise().norm().sum();
            std::cout<<(knots-knots_prev).colwise().norm().sum()/knots.colwise().norm().sum()<<std::endl;
            //std::cout << uh << std::endl;
            knots_prev = knots;
            std::cout << "length: " << lencomp()<<std::endl;
            
            q++;
        } while (q<5000 && it_err>1e-14);
        std::cout << knots << std::endl;
        std::cout << "finished\n";

        for(int i=0; i<M+1;++i) knots(0,i) = (M_PI*i)/M-std::sin((M_PI*i)/M);
        for(int i=0; i<M+1;++i) knots(1,i) = std::cos((M_PI*i)/M)-1.;
        std::cout << "length cycloid: " << lencomp()<<std::endl;
        return knots;
    };

}  // namespace Brachistochrone
