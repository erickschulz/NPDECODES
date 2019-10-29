#include <gtest/gtest.h>
#include "../local_assembler.h"

TEST(ElementMatrixTest, TestMass) {
  TriGeo_t elem;
  elem << 0., 1., 0., 0., 0., 1.;
  Eigen::MatrixXd element_matrix = ElementMatrix_Mass_LFE(elem);
  Eigen::MatrixXd correct_element_matrix(3, 3);
  correct_element_matrix << 0.0833333, 0.0416667, 0.0416667, 0.0416667,
      0.0833333, 0.0416667, 0.0416667, 0.0416667, 0.0833333;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(element_matrix(i, j), correct_element_matrix(i, j), 0.0001);
    }
  }
}

TEST(ElementMatrixTest, TestStiffness) {
  TriGeo_t elem;
  elem << 0., 1., 0., 0., 0., 1.;
  Eigen::MatrixXd element_matrix = ElementMatrix_Lapl_LFE(elem);
  Eigen::MatrixXd correct_element_matrix(3, 3);
  std::cout << element_matrix;
  correct_element_matrix << 1, -0.5, -0.5, -0.5, 0.5, 0, -0.5, 0, 0.5;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(element_matrix(i, j), correct_element_matrix(i, j), 0.0001);
    }
  }
}

TEST(ElementMatrixTest, TestLoad) {
  TriGeo_t elem;
  elem << 0., 1., 0., 0., 0., 1.;
  auto f = [](const Eigen::Vector2d &x) { return 1.; };
  Eigen::VectorXd element_vector = localLoadLFE(elem, f);
  EXPECT_NEAR(element_vector.sum(), 0.5, 0.0001);
}
