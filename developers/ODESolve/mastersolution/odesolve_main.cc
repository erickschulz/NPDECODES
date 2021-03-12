#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "odesolve.h"


int main()
{
    auto f = [] (const Eigen::VectorXd &y) -> Eigen::VectorXd { return Eigen::VectorXd::Ones(1) + y*y; };
    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double T = 1.;
    auto y_ex = [] (double t) -> Eigen::VectorXd { Eigen::VectorXd y(1); y << tan(t); return y; };
    
    unsigned p = 1;
    
    auto Psi = [&f] (double h, const Eigen::VectorXd & y0) -> Eigen::VectorXd { return y0 + h*f(y0); };
    auto PsiTilde = [&Psi, &f, &p] (double h, const Eigen::VectorXd & y0) -> Eigen::VectorXd {  return ODESolve::psitilde(Psi, p, h, y0); };
    
/* SAM_LISTING_BEGIN_2 */
    //// Subproblem d
    std::cout << "*** SUBPROBLEM d:" << std::endl;
    
#if SOLUTION
    std::cout << "Error table for equidistant steps:" << std::endl;
    std::cout << "N" << "\t" << "Error" << std::endl;
    for(int N = 4; N < 4096; N=N<<1 ) {
        std::vector<Eigen::VectorXd> Y = ODESolve::odeintequi(PsiTilde, T, y0, N);
        double err = (Y.back() - y_ex(1)).norm();
        std::cout << N << "\t" << err << std::endl;
    }
#else // TEMPLATE
    // TODO: determine order of convergence of PsiTilde using odeintequi
#endif // TEMPLATE
/* SAM_LISTING_END_2 */
    
    //// Subproblem e
    std::cout << "*** SUBPROBLEM e:" << std::endl;
    
    double h0 = 1./100.;
    std::vector<Eigen::VectorXd> Y;
    double err = 0;

#if SOLUTION
    Y = ODESolve::odeintssctrl(Psi, T, y0, h0, p, 10e-4, 10e-4, 10e-5);
    err = (Y.back() - y_ex(1)).norm();
    std::cout << "Adaptive error control results:" << std::endl;
    std::cout << "Error" << "\t" << "No. Steps" << "\t" << "y(1)" << "\t" << "y(ex)" << std::endl;
    std::cout << err << "\t" << Y.size() << "\t" << Y.back() << "\t" << y_ex(1) << std::endl;
#else // TEMPLATE
    // TODO: uncomment following lines to test ODESolve::odeintssctrl
//    Y = ODESolve::odeintssctrl(Psi, T, y0, h0, p, 10e-4, 10e-4, 10e-5);
//    err = (Y.back() - y_ex(1)).norm();
//    std::cout << "Adaptive error control results:" << std::endl;
//    std::cout << "Error" << "\t\t" << "No. Steps" << "\t" << "y(1)" << "\t" << "y(ex)" << std::endl;
//    std::cout << err << "\t" << Y.size() << "\t" << Y.back() << "\t" << y_ex(1) << std::endl;
#endif // TEMPLATE

    return 0;
}
