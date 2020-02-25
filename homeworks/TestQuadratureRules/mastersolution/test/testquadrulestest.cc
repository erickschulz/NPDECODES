/**
 * @file
 * @brief UNITESTS for NPDE homework TestQuadratureRules
 * @author Erick Schulz, Liaowang Huang (refactoring)
 * @date 08/03/2019, , 22/02/2020 (refactoring)
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>

#include <lf/base/base.h>
#include <lf/quad/quad.h>

#include "../testquadraturerules.h"

namespace TestQuadratureRules::test {

using namespace TestQuadratureRules;

TEST(TestQuadratureRules, TestQuadratureTria) {
  std::cout << "*********************************************" << std::endl;
  std::cout << "NPDE homework TestQuadratureRules: UNIT TESTS" << std::endl;
  std::cout << "*********************************************" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << "Testing quadratures over reference triangle." << std::endl;
  std::cout << "--------------------------------------------" << std::endl;


  constexpr unsigned int max_test_order = 15;  // This CANNOT be changed

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
      ASSERT_TRUE(result_mat_expected[i][j] == result_mat_kTria[i][j])
        << "Fail!" << std::endl;
    }
  }
  std::cout << "Pass!" << std::endl;
}

TEST(TestQuadratureRules, TestQuadratureQuad) {
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << "Testing quadratures over the reference cube." << std::endl;
  std::cout << "--------------------------------------------" << std::endl;

  constexpr unsigned int max_test_order = 15;  // This CANNOT be changed

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
      ASSERT_TRUE(result_mat_expected[i][j] == result_mat_kQuad[i][j]) << "Fail!" << std::endl;
    }
  }
  std::cout << "Pass!" << std::endl;
}

TEST(TestQuadratureRules, calcQuadOrder) {
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "           Testing calcQuadOrder            " << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    constexpr unsigned int max_test_order = 14;
    
    // the expected results to be compared.
    int result_tria_expected[max_test_order] = {2, 2, 3, 5, 5, 6, 7, 8, 9, 10, 11, 12,
        13, 14};
    int result_quad_expected[max_test_order] = {2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12,
        14, 14};
    
    // test calcQuadorder for reference triangle and quad
    const auto ref_tria = lf::base::RefEl(lf::base::RefElType::kTria);
    const auto ref_quad = lf::base::RefEl(lf::base::RefElType::kQuad);
    lf::quad::QuadRule quad_rule;

    // Create the quadratures and proceed to testing
    for (int order_requested = 1; order_requested <= max_test_order;
         order_requested++) {
        quad_rule = lf::quad::make_QuadRule(ref_tria, order_requested - 1);

        ASSERT_TRUE(calcQuadOrder(quad_rule) == result_tria_expected[order_requested-1])
        << "Fail!";
    }
    
    for (int order_requested = 1; order_requested <= max_test_order;
         order_requested++) {
        quad_rule = lf::quad::make_QuadRule(ref_quad, order_requested - 1);

        ASSERT_TRUE(calcQuadOrder(quad_rule) == result_quad_expected[order_requested-1])
        << "Fail!" << std::endl;
    }
    std::cout << "Pass!" << std::endl;
}

}  // namespace TestQuadratureRules::test
