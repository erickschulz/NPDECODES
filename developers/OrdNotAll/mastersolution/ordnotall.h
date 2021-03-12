#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "rkintegrator.h"


namespace OrdNotAll {

/*!
 * \brief errors Approximates the order of convergence of a scheme.
 * The scheme is a Runge Kutta scheme defined by A and b when applied
 * to the first order system y'=f(y), y(0)=y0.
 * The function ouputs error of the solutions at the time T.
 * \tparam Function Type for r.h.s function f.
 * \param f The r.h.s function for the ODE.
 * \param T Final time.
 * \param y0 Initial data.
 * \param A Butcher matrix $A$.
 * \param b Buthcer vector $b$.
 */
/* SAM_LISTING_BEGIN_1 */
template <class Function>
void errors(const Function &f, double T,
            const Eigen::VectorXd &y0,
            const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {

    RKIntegrator<Eigen::VectorXd> rk(A, b);

    std::vector<double> error(15);
#if SOLUTION
    std::vector<double> order(14);

    double sum = 0;
    int count = 0;
    bool test = 1;

    std::vector<Eigen::VectorXd> y_exact = rk.solve(f, T, y0, std::pow(2,15));

    for(int k = 0; k < 15; k++) {
        int N = std::pow(2,k+1);
        std::vector<Eigen::VectorXd> y1 = rk.solve(f,T,y0,N);

        error[k] = (y1[N] - y_exact[std::pow(2,15)]).norm();

        std::cout << std::left << std::setfill(' ')
                  << std::setw(3)  << "N = "
                  << std::setw(7) << N
                  << std::setw(8) << "Error = "
                  << std::setw(13) << error[k];

        if ( error[k] < y0.size() * 5e-14) {
            test = 0;
        }
        if ( k>0 && test ) {
            order[k-1]=std::log2( error[k-1] / error[k] );
            std::cout << std::left << std::setfill(' ')
                      << std::setw(10)
                      << "Approximated order = " << order[k-1]
                      << std::endl;
            sum += order[k-1];
            count = k;
        }
        else std::cout << std::endl;
    }
    std::cout << "Average approximated order = "
              << sum / count
              << std::endl << std::endl;
#else // TEMPLATE
    // TODO: output error and order of the method
#endif
}
/* SAM_LISTING_END_1 */

}  // namespace OrdNotAll
