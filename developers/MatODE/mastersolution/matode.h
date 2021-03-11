#include <limits>

#include <Eigen/Dense>

#include "ode45.h"


namespace MatODE {

//! \brief Solve matrix IVP Y' = -(Y-Y')*Y using ode45 up to time T
//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
//! \return Matrix of solution of IVP at t = T
/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd matode(const Eigen::MatrixXd & Y0, double T) {
#if SOLUTION
    auto F = [] (const Eigen::MatrixXd & M) {
        return -(M  - M.transpose())*M;
    };
    ode45<Eigen::MatrixXd> O(F);

    // Set tolerances
    O.options.atol = 10e-10;
    O.options.rtol = 10e-8;

    // Return only matrix at $T$, (solution is vector
    // of pairs $(y(t_k), t_k)$ for each step k
    return O.solve(Y0, T).back().first;
#else // TEMPLATE
    // TODO: solve matrix ODE with ode45 class
    return Y0;
#endif // TEMPLATE
}
/* SAM_LISTING_END_1 */

//! \brief Find if invariant is preserved after evolution with matode
//! \param[in] Y0 Initial data Y(0) (as matrix)
//! \param[in] T final time of simulation
//! \return true if invariant was preserved (up to round-off),
//! i.e. if norm was less than 10*eps
/* SAM_LISTING_BEGIN_2 */
bool checkinvariant(const Eigen::MatrixXd & M, double T) {
#if SOLUTION
    Eigen::MatrixXd N(3,3);

    N = matode(M, T);

    if( (N.transpose()*N-M.transpose()*M).norm() <
        10 * std::numeric_limits<double>::epsilon()* M.norm()) {
        return true;
    } else {
        return false;
    }
#else // TEMPLATE
    // TODO: test wether matrix satisfy invariant.
    return false;
#endif // TEMPLATE
}
/* SAM_LISTING_END_2 */

//! \brief Implement ONE step of explicit Euler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd expeulstep(const Eigen::MatrixXd & A, const Eigen::MatrixXd & Y0, double h) {
#if SOLUTION
    return Y0 + h*A*Y0;
#else // TEMPLATE
    // TODO: implement one explicit Euler step
    return Y0;
#endif // TEMPLATE
}
/* SAM_LISTING_END_3 */

//! \brief Implement ONE step of implicit Euler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
/* SAM_LISTING_BEGIN_4 */
Eigen::MatrixXd impeulstep(const Eigen::MatrixXd & A, const Eigen::MatrixXd & Y0, double h) {
#if SOLUTION
    return (Eigen::MatrixXd::Identity(3,3) - h*A).partialPivLu().solve(Y0);
#else // TEMPLATE
    // TODO: implement one implicit Euler step
    return Y0;
#endif // TEMPLATE
}
/* SAM_LISTING_END_4 */

//! \brief Implement ONE step of implicit midpoint ruler applied to Y0, of ODE Y' = A*Y
//! \param[in] A matrix A of the ODE
//! \param[in] Y0 Initial state
//! \param[in] h step size
//! \return next step
/* SAM_LISTING_BEGIN_5 */
Eigen::MatrixXd impstep(const Eigen::MatrixXd & A, const Eigen::MatrixXd & Y0, double h) {
#if SOLUTION
    return (Eigen::MatrixXd::Identity(3,3) - h*0.5*A)
            .partialPivLu()
            .solve(Y0+h*0.5*A*Y0);
#else // TEMPLATE
    // TODO: implement one implicit midpoint step.
    return Y0;
#endif // TEMPLATE
}
/* SAM_LISTING_END_5 */

}  // namespace MatODE
