#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/QR>

#include "matode.h"


int main() {
    /* SAM_LISTING_BEGIN_6 */
    double T = 1;
    unsigned int n = 3;

    Eigen::MatrixXd M(n,n);
    M << 8,1,6,3,5,7,4,9,2;

    std::cout << "SUBTASK c)" << std::endl;
    // Test preservation of orthogonality

    // Build Q
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(M.rows(), M.cols());
    qr.compute(M);
    Eigen::MatrixXd Q = qr.householderQ();

    // Build A
    Eigen::MatrixXd A(n,n);
    A << 0, 1, 1, -1, 0, 1, -1, -1, 0;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n,n);

#if SOLUTION
    // Norm of Y'Y-I for 20 steps
    Eigen::MatrixXd Mexpeul = Q, Mimpeul = Q, Mimp = Q;
    double h = 0.01;
    std::vector<int> sep = {8,15,15,15};
    std::cout << "Evolution of norm(Y_k'*Y_k - I) for three methods:" << std::endl;
    std::cout   << std::setw(sep[0]) << "step"
                << std::setw(sep[1]) << "exp. Eul"
                << std::setw(sep[2]) << "imp. Eul"
                << std::setw(sep[3]) << "IMP"
                << std::endl;
    std::cout   << std::setw(sep[0]) << "-1"
                << std::setw(sep[1]) << (Mexpeul.transpose()*Mexpeul - I).norm()
                << std::setw(sep[2]) << (Mimpeul.transpose()*Mimpeul - I).norm()
                << std::setw(sep[3]) << (Mimp.transpose()*Mimp - I).norm()
                << std::endl;
    for(unsigned int j = 0; j < 20; ++j) {
        Mexpeul = MatODE::expeulstep(A, Mexpeul, h);
        Mimpeul = MatODE::impeulstep(A, Mimpeul, h);
        Mimp = MatODE::impstep(A, Mimp, h);

        std::cout   << std::setw(sep[0]) << j
                    << std::setw(sep[1]) << (Mexpeul.transpose()*Mexpeul - I).norm()
                    << std::setw(sep[2]) << (Mimpeul.transpose()*Mimpeul - I).norm()
                    << std::setw(sep[3]) << (Mimp.transpose()*Mimp - I).norm()
                    << std::endl;
    }
#else // TEMPLATE
    // TODO
#endif // TEMPLATE
    /* SAM_LISTING_END_6 */

    /* SAM_LISTING_BEGIN_7 */
    std::cout << "Test implementation of ode45" << std::endl;
#if SOLUTION
    std::cout << "M = " << std::endl
              << M << std::endl;
    Eigen::MatrixXd  N = MatODE::matode(M, T);
    std::cout << "N = " << std::endl
              << N << std::endl;
#else // TEMPLATE
    // TODO
#endif // TEMPLATE
    /* SAM_LISTING_END_7 */

    /* SAM_LISTING_BEGIN_8 */
    std::cout << "Test whether invariant was preserved or not" << std::endl;
#if SOLUTION
    bool is_invariant = MatODE::checkinvariant(N, T);

    if( is_invariant ) {
        std::cout << "Invariant was preserved."
                  << std::endl;
    } else {
        std::cout << "Invariant was NOT preserved."
                  << std::endl;
    }
#else // TEMPLATE
    // TODO
#endif // TEMPLATE
    /* SAM_LISTING_END_8 */
}
