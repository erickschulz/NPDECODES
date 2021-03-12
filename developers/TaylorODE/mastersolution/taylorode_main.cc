#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "taylorode.h"

//! \file taylorprey.cpp Problem 3, solving prey/predator model with taylor expansion method

int main() {
/* SAM_LISTING_BEGIN_0 */
    // Dimension of state space
    unsigned int d = 2;
    
    // Final time for model
    double T = 10.;
    
    // Initial value for model
    Eigen::VectorXd y0(d);
    y0 << 100, 5;
    
    // Array of number of steps (for convergence study)
    std::vector<unsigned int> N = {128, 256, 512, 1024, 2048, 4096, 8192, 16384};
    
    // Exact value $y(10)$ at final time $T = 10$ (approximated)
    Eigen::VectorXd yex(d);
    yex << 0.319465882659820, 9.730809352326228;
    
#if SOLUTION
    // Structure for the evaluation of $f$, $df$ and $d2f$ (Hessian bilinear form)
    struct odefun {
        std::function<Eigen::VectorXd(const Eigen::VectorXd &)> f = [this] (const Eigen::VectorXd & y) {
            Eigen::VectorXd temp = y;
            temp(0) *= (alpha1 - beta1*y(1));
            temp(1) *= (beta2*y(0) - alpha2);
            return temp;
        };
        std::function<Eigen::VectorXd(const Eigen::VectorXd &, const Eigen::VectorXd &)> df = [this] (const Eigen::VectorXd & y, const Eigen::VectorXd & z) {
            Eigen::MatrixXd Df(2,2);
            Df << alpha1-beta1*y(1), -beta1*y(0),  beta2*y(1), -alpha2+beta2*y(0);
            Eigen::VectorXd temp = Df*z;
            return temp;
        };
        std::function<Eigen::VectorXd(const Eigen::VectorXd &, const Eigen::VectorXd &)> d2f = [this] (const Eigen::VectorXd & y, const Eigen::VectorXd & z) {
            Eigen::MatrixXd H1(2,2), H2(2,2);
            H1 << 0, -beta1, -beta1, 0;
            H2 << 0, beta2, beta2, 0;
            Eigen::VectorXd temp(2);
            temp(0) = z.transpose()*H1*z;
            temp(1) = z.transpose()*H2*z;
            return temp;
        };
        const double alpha1 = 3;
        const double alpha2 = 2;
        const double beta1 = 0.1;
        const double beta2 = 0.1;
    } F;
#else // TEMPLATE
    // Coefficients and handle for prey/predator model
    double alpha1 = 3.;
    double alpha2 = 2.;
    double beta1 = 0.1;
    double beta2 = 0.1;
    // TODO: implement a structure for the evaluation of $f$, $df$ and $d2f$ (Hessian bilinear form)
#endif // TEMPLATE
    
    // Constructor
    TaylorODE::TaylorIntegrator<Eigen::VectorXd> tint;
    
    // Start convergence study
    std::cout  << std::setw(15) << "N"  << std::setw(15) << "error" << std::setw(15) << "rate" << std::endl;
    #if SOLUTION
    double errold;
    for(unsigned int i = 0; i < N.size(); ++i) {
        auto res = tint.solve(F, T, y0, N[i]);
        double err = (res.back() - yex).norm();
        std::cout  << std::setw(15) << N[i] << std::setw(15) << err;
        if(i > 0) {
            std::cout << std::setw(15) << std::log2(errold / err);
        }
        errold = err;
        std::cout  << std::endl;
    }
    #else // TEMPLATE
        // TODO: tabulate error for each $N$
    #endif // TEMPLATE
/* SAM_LISTING_END_0 */

    return 0;
}
