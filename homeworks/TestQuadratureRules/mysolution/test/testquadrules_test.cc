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
      EXPECT_TRUE(result_mat_expected[i][j] == result_mat_kTria[i][j]);
    }
  }
}

TEST(TestQuadratureRules, TestQuadratureQuad) {
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
      EXPECT_TRUE(result_mat_expected[i][j] == result_mat_kQuad[i][j]);
    }
  }
}

TEST(TestQuadratureRules, calcQuadOrder) {
  constexpr unsigned int max_test_order = 30;

  // the expected results to be compared.
  int result_tria_expected[max_test_order] = {
      2,  2,  3,  5,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
      16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
  int result_quad_expected[max_test_order] = {
      2,  2,  4,  4,  6,  6,  8,  8,  10, 10, 12, 12, 14, 14, 16,
      16, 18, 18, 20, 20, 22, 22, 25, 25, 28, 28, 32, 32, 35, 35};

  // test calcQuadorder for reference triangle and quad
  const auto ref_tria = lf::base::RefEl(lf::base::RefElType::kTria);
  const auto ref_quad = lf::base::RefEl(lf::base::RefElType::kQuad);
  lf::quad::QuadRule quad_rule;

  // Create the quadratures and proceed to testing
  for (int order_requested = 1; order_requested <= max_test_order;
       order_requested++) {
    quad_rule = lf::quad::make_QuadRule(ref_tria, order_requested - 1);

    EXPECT_TRUE(calcQuadOrder(quad_rule) ==
                result_tria_expected[order_requested - 1]);
  }

  for (int order_requested = 1; order_requested <= max_test_order;
       order_requested++) {
    quad_rule = lf::quad::make_QuadRule(ref_quad, order_requested - 1);

    EXPECT_TRUE(calcQuadOrder(quad_rule) ==
                result_quad_expected[order_requested - 1]);
  }
}

}  // namespace TestQuadratureRules::test
