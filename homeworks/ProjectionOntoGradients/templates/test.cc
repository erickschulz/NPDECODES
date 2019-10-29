#include <gtest/gtest.h>
#include "gradprojection.h"

TEST(GradProjection, div_free_test) {
  EXPECT_EQ(0, 42) << "You still need to implement this test. Remove this line "
                      "once you start implementing ";
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION
}

TEST(GradProjection, exact_sol_test) {
  // mesh builder in a world of dimension 2
  lf::mesh::hybrid2d::TPTriagMeshBuilder my_builder(
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2));

  // define the test mesh
  my_builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNoXCells(10)
      .setNoYCells(10);

  auto mesh_p = my_builder.Build();

  EXPECT_EQ(0, 42) << "You still need to implement this test. Remove this line "
                      "once you start implementing ";

  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION
}
