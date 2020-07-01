/**
 * @file
 * @brief NPDE homework TestQuadratureRules
 * @author Erick Schulz, Liaowang Huang (refactoring)
 * @date 08/03/2019, 22/02/2020 (refactoring)
 * @copyright Developed at ETH Zurich
 */

#include <iostream>

#include <lf/base/base.h>
#include <lf/quad/quad.h>

#include "testquadraturerules.h"

using namespace TestQuadratureRules;

int main() {
  std::cout << " " << std::endl;
  std::cout << "********************************************" << std::endl;
  std::cout << " NPDE Homework problem: TestQuadratureRules " << std::endl;
  std::cout << "********************************************" << std::endl;

  unsigned int max_order_tested = 30;
  unsigned int order;

  std::cout << " " << std::endl;
  std::cout << "Testing orders of LehrFEM++ quadratures." << std::endl;
  std::cout << " " << std::endl;
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "|                RESULTS                  |" << std::endl;

  std::vector<int> max_orders_kTria;
  std::vector<int> max_orders_kQuad;

  for (order = 1; order <= max_order_tested; order++) {
    // Quadrature for reference triangle
    const auto ref_triangle = lf::base::RefEl(lf::base::RefElType::kTria);
    const auto quad_rule_kTria =
        lf::quad::make_QuadRule(ref_triangle, order - 1);

    // Quadrature for reference cube
    const auto ref_cube = lf::base::RefEl(lf::base::RefElType::kQuad);
    const auto quad_rule_kQuad = lf::quad::make_QuadRule(ref_cube, order - 1);

    // Testing for maximal order
    unsigned int cur_maximal_order_kTria = calcQuadOrder(quad_rule_kTria);
    max_orders_kTria.push_back(cur_maximal_order_kTria);
    unsigned int cur_maximal_order_kQuad = calcQuadOrder(quad_rule_kQuad);
    max_orders_kQuad.push_back(cur_maximal_order_kQuad);
  }

  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "|       |          MAXIMAL ORDER          |" << std::endl;
  std::cout << "| ORDER |---------------------------------|" << std::endl;
  std::cout << "|       |     kTria     |       kQuad     |" << std::endl;
  std::cout << "-------------------------------------------" << std::endl;
  for (order = 1; order <= max_order_tested; order++) {
    if (order < 9) {
      std::cout << "|     " << order << " |"
                << "       " << max_orders_kTria[order - 1] << "    "
                << "   |      " << max_orders_kQuad[order - 1] << "          |"
                << std::endl;
    } else if (order == 9) {
      std::cout << "|     " << order << " |"
                << "       " << max_orders_kTria[order - 1] << "    "
                << "   |      " << max_orders_kQuad[order - 1] << "         |"
                << std::endl;
    } else {
      std::cout << "|    " << order << " |"
                << "       " << max_orders_kTria[order - 1] << "    "
                << "  |      " << max_orders_kQuad[order - 1] << "         |"
                << std::endl;
    }
  }
  std::cout << "-------------------------------------------" << std::endl;

  return 0;
}
