/**
 * @file
 * @brief NPDE homework TestQuadratureRules
 * @author Erick Schulz
 * @date 08/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "test_quad_rules.h"

namespace TestQuadratureRules {

using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;

bool testQuadOrderTria(const lf::quad::QuadRule &quad_rule,
                       unsigned int order) {
  bool order_isExact = true;  // return variable
  /* SOLUTION_BEGIN */
  
  /* SOLUTION_END */
  return order_isExact;
}

bool testQuadOrderQuad(const lf::quad::QuadRule &quad_rule,
                       unsigned int order) {
  bool order_isExact = true;  // return variable
  /* SOLUTION_BEGIN */
  
  /* SOLUTION_END */
  return order_isExact;
}

unsigned int calcQuadOrder(const lf::quad::QuadRule &quad_rule) {
  unsigned int maximal_order = quad_rule.Order();

  /* SOLUTION_BEGIN */

  /* SOLUTION_END */
  return maximal_order;
}

}  // namespace TestQuadratureRules
