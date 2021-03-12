#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include "systemode.h"


int main() {
/* SAM_LISTING_BEGIN_0 */
#if SOLUTION
    // Construct data for RK order 4
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(4,4);
    A(1,0) = .5;
    A(2,1) = .5;
    A(3,2) = 1;
    Eigen::VectorXd b(4);
    b << 1./6, 1./3, 1./3, 1./6;
    
    // Construct data for the IVP
    double T = 1;
    int n = 5;
    Eigen::VectorXd y0(2*n);
    for(int i = 0; i < n; ++i) {
        y0(i)=(i+1.)/n;
        y0(i+n)=-1;
    }

    auto f = [n] (Eigen::VectorXd y) {
        Eigen::VectorXd fy(2*n);
        
        Eigen::VectorXd g(n);
        g(0) = y(0)*(y(1)+y(0));
        g(n-1) = y(n-1)*(y(n-1)+y(n-2));
        for(int i = 1; i < n-1; ++i) {
            g(i) = y(i)*(y(i-1)+y(i+1));
        }
        
        Eigen::SparseMatrix<double> C(n,n);
        C.reserve(3);
        for(int i = 0; i < n; ++i) {
            C.insert(i,i) = 2;
            if(i < n-1) C.insert(i,i+1) = -1;
            if(i >= 1)  C.insert(i,i-1) = -1;
        }
        C.makeCompressed();
//      fy.head(n) = y.head(n);
        fy.head(n) = y.tail(n);
        
        Eigen::SparseLU< Eigen::SparseMatrix<double> > solver;
        solver.analyzePattern(C);
        solver.compute(C);
        fy.tail(n) = solver.solve(g);
        return fy;
    };

    errors(f, T, y0, A, b);
#else // TEMPLATE
    // TODO: compute $f$, $T$, $y_0$, $\VA$, $\Vb$ to run "errors(f, T, y0, A, b);"
#endif // TEMPLATE
/* SAM_LISTING_END_0 */

    return 0;
}
