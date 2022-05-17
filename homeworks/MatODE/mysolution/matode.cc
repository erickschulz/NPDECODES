/**
 * @file matode.cc
 * @brief NPDE homework MatODE code
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>

namespace MatODE {

/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd eeulstep(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Y0,
                         double h) {
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(Y0.rows(), Y0.cols());
  //====================
  // Your code goes here
  //====================
  return res;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::MatrixXd ieulstep(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Y0,
                         double h) {
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(Y0.rows(), Y0.cols());
  //====================
  // Your code goes here
  //====================
  return res;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
Eigen::MatrixXd impstep(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Y0,
                        double h) {
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(Y0.rows(), Y0.cols());
  //====================
  // Your code goes here
  //====================
  return res;
}
/* SAM_LISTING_END_5 */

}  // namespace MatODE
