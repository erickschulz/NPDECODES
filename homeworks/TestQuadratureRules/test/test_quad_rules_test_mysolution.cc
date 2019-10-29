/**
 * @file
 * @brief UNITESTS for NPDE homework TestQuadratureRules
 * @author Erick Schulz
 * @date 08/03/2019
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>

#include "../mysolution/test_quad_rules.h"

namespace TestQuadratureRules::test {

using namespace TestQuadratureRules;

TEST(TestQuadratureRules, TestQuadratureTria) {
  std::cout << "" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "*********************************************" << std::endl;
  std::cout << "NPDE homework TestQuadratureRules: UNIT TESTS" << std::endl;
  std::cout << "*********************************************" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << "Testing quadratures over reference triangle." << std::endl;
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << "Displaying obtained ouput:"
            << "\n"
            << std::endl;

  constexpr unsigned int max_test_order = 15;  // This CANNOT be changed

  // Some formatting for the ouput table
  for (int k = 1; k <= max_test_order; k++) {
    if (k < 10) {
      std::cout << k << "     ";
    } else {
      std::cout << k << "    ";
    }
  }

  // Result boolean double array whose columns indices specify tested order
  // and rows indices the expected order of the quadrature
  bool result_mat_kTria[max_test_order][max_test_order];

  // Quadrature for reference triangle
  const auto ref_triangle = lf::base::RefEl(lf::base::RefElType::kTria);
  lf::quad::QuadRule quad_rule;

  // Create the quadratures and proceed to testing
  for (int order_requested = 1; order_requested <= max_test_order;
       order_requested++) {
    quad_rule = lf::quad::make_QuadRule(ref_triangle, order_requested - 1);

    for (int order_tested = 1; order_tested <= max_test_order; order_tested++) {
      result_mat_kTria[order_tested - 1][order_requested - 1] =
          testQuadOrderTria(quad_rule, order_tested);
    }
  }

  // Display the obtained results to the terminal
  std::cout << " " << std::endl;
  for (int order_tested = 0; order_tested < max_test_order; order_tested++) {
    for (int order_requested = 0; order_requested < max_test_order;
         order_requested++) {
      if (result_mat_kTria[order_tested][order_requested]) {
        std::cout << std::boolalpha
                  << result_mat_kTria[order_tested][order_requested];
        std::cout << " ";
        std::cout << " ";
      } else {
        std::cout << std::boolalpha
                  << result_mat_kTria[order_tested][order_requested];
        std::cout << " ";
      }
    }
    printf("\n");
  }
  printf("\n");
  printf("\n");

  // Compare the results to the expected results
  bool result_mat_expected[max_test_order][max_test_order] = {
      {true, true, true, true, true, true, true, true, true, true, true, true,
       true, true, true},
      {true, true, true, true, true, true, true, true, true, true, true, true,
       true, true, true},
      {false, false, true, true, true, true, true, true, true, true, true, true,
       true, true, true},
      {false, false, false, true, true, true, true, true, true, true, true,
       true, true, true, true},
      {false, false, false, true, true, true, true, true, true, true, true,
       true, true, true, true},
      {false, false, false, false, false, true, true, true, true, true, true,
       true, true, true, true},
      {false, false, false, false, false, false, true, true, true, true, true,
       true, true, true, true},
      {false, false, false, false, false, false, false, true, true, true, true,
       true, true, true, true},
      {false, false, false, false, false, false, false, false, true, true, true,
       true, true, true, true},
      {false, false, false, false, false, false, false, false, false, true,
       true, true, true, true, true},
      {false, false, false, false, false, false, false, false, false, false,
       true, true, true, true, true},
      {false, false, false, false, false, false, false, false, false, false,
       false, true, true, true, true},
      {false, false, false, false, false, false, false, false, false, false,
       false, false, true, true, true},
      {false, false, false, false, false, false, false, false, false, false,
       false, false, false, true, true},
      {false, false, false, false, false, false, false, false, false, false,
       false, false, false, false, true},
  };

  for (int i = 0; i < max_test_order; i++) {
    for (int j = 0; j < max_test_order; j++) {
      ASSERT_TRUE(result_mat_expected[i][j] == result_mat_kTria[i][j]);
    }
  }
}

TEST(TestQuadratureRules, TestQuadratureQuad) {
  std::cout << "" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << "Testing quadratures over the reference cube." << std::endl;
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << "Displaying obtained ouput:"
            << "\n"
            << std::endl;

  constexpr unsigned int max_test_order = 15;  // This CANNOT be changed

  // Some formatting for the ouput table
  for (int k = 1; k <= max_test_order; k++) {
    if (k < 10) {
      std::cout << k << "     ";
    } else {
      std::cout << k << "    ";
    }
  }

  // Result boolean double array whose columns indices specify tested order
  // and rows indices the expected order of the quadrature
  bool result_mat_kQuad[max_test_order][max_test_order];

  // Quadrature for reference cube
  const auto ref_cube = lf::base::RefEl(lf::base::RefElType::kQuad);
  lf::quad::QuadRule quad_rule;

  // Create the quadratures and proceed to testing
  for (int order_requested = 1; order_requested <= max_test_order;
       order_requested++) {
    quad_rule = lf::quad::make_QuadRule(ref_cube, order_requested - 1);

    for (int order_tested = 1; order_tested <= max_test_order; order_tested++) {
      result_mat_kQuad[order_tested - 1][order_requested - 1] =
          testQuadOrderQuad(quad_rule, order_tested);
    }
  }

  // Display the obtained results to the terminal
  std::cout << " " << std::endl;
  for (int order_tested = 0; order_tested < max_test_order; order_tested++) {
    for (int order_requested = 0; order_requested < max_test_order;
         order_requested++) {
      if (result_mat_kQuad[order_tested][order_requested]) {
        std::cout << std::boolalpha
                  << result_mat_kQuad[order_tested][order_requested];
        std::cout << " ";
        std::cout << " ";
      } else {
        std::cout << std::boolalpha
                  << result_mat_kQuad[order_tested][order_requested];
        std::cout << " ";
      }
    }
    printf("\n");
  }
  printf("\n");
  printf("\n");

  // Compare the results to the expected results
  bool result_mat_expected[max_test_order][max_test_order] = {
      {true, true, true, true, true, true, true, true, true, true, true, true,
       true, true, true},
      {true, true, true, true, true, true, true, true, true, true, true, true,
       true, true, true},
      {false, false, true, true, true, true, true, true, true, true, true, true,
       true, true, true},
      {false, false, true, true, true, true, true, true, true, true, true, true,
       true, true, true},
      {false, false, false, false, true, true, true, true, true, true, true,
       true, true, true, true},
      {false, false, false, false, true, true, true, true, true, true, true,
       true, true, true, true},
      {false, false, false, false, false, false, true, true, true, true, true,
       true, true, true, true},
      {false, false, false, false, false, false, true, true, true, true, true,
       true, true, true, true},
      {false, false, false, false, false, false, false, false, true, true, true,
       true, true, true, true},
      {false, false, false, false, false, false, false, false, true, true, true,
       true, true, true, true},
      {false, false, false, false, false, false, false, false, false, false,
       true, true, true, true, true},
      {false, false, false, false, false, false, false, false, false, false,
       true, true, true, true, true},
      {false, false, false, false, false, false, false, false, false, false,
       false, false, true, true, true},
      {false, false, false, false, false, false, false, false, false, false,
       false, false, true, true, true},
      {false, false, false, false, false, false, false, false, false, false,
       false, false, false, false, true}};

  for (int i = 0; i < max_test_order; i++) {
    for (int j = 0; j < max_test_order; j++) {
      ASSERT_TRUE(result_mat_expected[i][j] == result_mat_kQuad[i][j]);
    }
  }
}

}  // namespace TestQuadratureRules::test
