#include <gtest/gtest.h>
#include "gradprojection.h"

/* SAM_LISTING_BEGIN_1 */
TEST(GradProjection, div_free_test) {
  /* BEGIN_SOLUTION */

  // Initialize all objects we need for our test case
  const auto f = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x(1), x(0));
  };
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  // DofHandler for the p.w. linear Lagrangian finite element method
  lf::assemble::UniformFEDofHandler dofh(mesh_p,
                                         {{lf::base::RefEl::kPoint(), 1},
                                          {lf::base::RefEl::kSegment(), 0},
                                          {lf::base::RefEl::kTria(), 0}});
  // Compute solution
  const Eigen::VectorXd sol_vec =
      ProjectionOntoGradients::projectOntoGradients(dofh, f);

  // As stated in the exercise, we would expect the solution to be zero.
  // Hence we check every entry if its (numerically) zero.
  // The GoogleTest framework provides a function we may use:
  /*   EXPECT_NEAR(value 1, value 2, max. difference) */
  const double eps = 1e-15;
  for (std::size_t i = 0; i < sol_vec.size(); ++i) {
    EXPECT_NEAR(sol_vec[i], 0.0, eps);
    // Try testing for equality, you'll see it will fail miserably!
    /* EXPECT_EQ(sol_vec[i], 0.0); */
  }
  /* END_SOLUTION */
}
/* SAM_LISTING_END_1 */

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

  /* BEGIN_SOLUTION */

  // Here we assign one degree of freedom to each node of our mesh
  // We assign 0 dof to edges and triangles
  // This corresponds to using the 'tent' basis functions of each node
  lf::assemble::UniformFEDofHandler dofh(mesh_p,
                                         {{lf::base::RefEl::kPoint(), 1},
                                          {lf::base::RefEl::kSegment(), 0},
                                          {lf::base::RefEl::kTria(), 0}});

  // define the function f
  const auto f = [&dofh](Eigen::Vector2d x) -> Eigen::Vector2d {
    int triang_idx = 9;
    // first determine of which triangle the coordinate is part of
    if (x(0) >= 0.0 && x(0) <= 0.5 && x(1) >= 0.0 && x(1) <= 0.5) {
      if (x(0) < x(1))
        triang_idx = 0;
      else
        triang_idx = 1;
    } else if (x(0) >= 0.0 && x(0) <= 0.5 && x(1) >= 0.5 && x(1) <= 1.0) {
      if (x(0) + 0.5 < x(1))
        triang_idx = 2;
      else
        triang_idx = 3;
    } else if (x(0) >= 0.5 && x(0) <= 1.0 && x(1) >= 0.0 && x(1) <= 0.5) {
      if (x(0) < 0.5 + x(1))
        triang_idx = 4;
      else
        triang_idx = 5;
    } else if (x(0) >= 0.5 && x(0) <= 1.0 && x(1) >= 0.5 && x(1) <= 1.0) {
      if (x(0) < x(1))
        triang_idx = 6;
      else
        triang_idx = 7;
    }
    // then return the function value according to the definition in the
    // exercise
    switch (triang_idx) {
      case 0:
        return Eigen::Vector2d(2, 0);
        break;
      case 1:
        return Eigen::Vector2d(0, 2);
        break;
      case 3:
        return Eigen::Vector2d(2, -2);
        break;
      case 4:
        return Eigen::Vector2d(-2, 2);
        break;
      case 6:
        return Eigen::Vector2d(0, -2);
        break;
      case 7:
        return Eigen::Vector2d(-2, 0);
        break;
      default:
        return Eigen::Vector2d(0, 0);
        break;
    }
  };

  // call the function you have written to test it
  const Eigen::VectorXd sol_vec =
      ProjectionOntoGradients::projectOntoGradients(dofh, f);

  // this is the residual we will accept as still correct
  const double eps = 1e-15;

  // hard-coded solution vector
  EXPECT_NEAR(sol_vec[0], 0, eps);
  EXPECT_NEAR(sol_vec[1], 0, eps);
  EXPECT_NEAR(sol_vec[2], 0, eps);
  EXPECT_NEAR(sol_vec[3], 0, eps);
  EXPECT_NEAR(sol_vec[4], 0, eps);
  EXPECT_NEAR(sol_vec[5], 0, eps);
  EXPECT_NEAR(sol_vec[6], 0, eps);
  EXPECT_NEAR(sol_vec[7], 0, eps);
  EXPECT_NEAR(sol_vec[8], 0, eps);
  EXPECT_NEAR(sol_vec[9], 0, eps);
  EXPECT_NEAR(sol_vec[10], 0, eps);
  EXPECT_NEAR(sol_vec[11], 0, eps);
  EXPECT_NEAR(sol_vec[12], 0.2, eps);
  EXPECT_NEAR(sol_vec[13], 0.2, eps);
  EXPECT_NEAR(sol_vec[14], 0.2, eps);
  EXPECT_NEAR(sol_vec[15], 0.2, eps);
  EXPECT_NEAR(sol_vec[16], 0.2, eps);
  EXPECT_NEAR(sol_vec[17], 0, eps);
  EXPECT_NEAR(sol_vec[18], 0, eps);
  EXPECT_NEAR(sol_vec[19], 0, eps);
  EXPECT_NEAR(sol_vec[20], 0, eps);
  EXPECT_NEAR(sol_vec[21], 0, eps);
  EXPECT_NEAR(sol_vec[22], 0, eps);
  EXPECT_NEAR(sol_vec[23], 0.2, eps);
  EXPECT_NEAR(sol_vec[24], 0.4, eps);
  EXPECT_NEAR(sol_vec[25], 0.4, eps);
  EXPECT_NEAR(sol_vec[26], 0.4, eps);
  EXPECT_NEAR(sol_vec[27], 0.4, eps);
  EXPECT_NEAR(sol_vec[28], 0.2, eps);
  EXPECT_NEAR(sol_vec[29], 0, eps);
  EXPECT_NEAR(sol_vec[30], 0, eps);
  EXPECT_NEAR(sol_vec[31], 0, eps);
  EXPECT_NEAR(sol_vec[32], 0, eps);
  EXPECT_NEAR(sol_vec[33], 0, eps);
  EXPECT_NEAR(sol_vec[34], 0.2, eps);
  EXPECT_NEAR(sol_vec[35], 0.4, eps);
  EXPECT_NEAR(sol_vec[36], 0.6, eps);
  EXPECT_NEAR(sol_vec[37], 0.6, eps);
  EXPECT_NEAR(sol_vec[38], 0.6, eps);
  EXPECT_NEAR(sol_vec[39], 0.4, eps);
  EXPECT_NEAR(sol_vec[40], 0.2, eps);
  EXPECT_NEAR(sol_vec[41], 0, eps);
  EXPECT_NEAR(sol_vec[42], 0, eps);
  EXPECT_NEAR(sol_vec[43], 0, eps);
  EXPECT_NEAR(sol_vec[44], 0, eps);
  EXPECT_NEAR(sol_vec[45], 0.2, eps);
  EXPECT_NEAR(sol_vec[46], 0.4, eps);
  EXPECT_NEAR(sol_vec[47], 0.6, eps);
  EXPECT_NEAR(sol_vec[48], 0.8, eps);
  EXPECT_NEAR(sol_vec[49], 0.8, eps);
  EXPECT_NEAR(sol_vec[50], 0.6, eps);
  EXPECT_NEAR(sol_vec[51], 0.4, eps);
  EXPECT_NEAR(sol_vec[52], 0.2, eps);
  EXPECT_NEAR(sol_vec[53], 0, eps);
  EXPECT_NEAR(sol_vec[54], 0, eps);
  EXPECT_NEAR(sol_vec[55], 0, eps);
  EXPECT_NEAR(sol_vec[56], 0.2, eps);
  EXPECT_NEAR(sol_vec[57], 0.4, eps);
  EXPECT_NEAR(sol_vec[58], 0.6, eps);
  EXPECT_NEAR(sol_vec[59], 0.8, eps);
  EXPECT_NEAR(sol_vec[60], 1, eps);
  EXPECT_NEAR(sol_vec[61], 0.8, eps);
  EXPECT_NEAR(sol_vec[62], 0.6, eps);
  EXPECT_NEAR(sol_vec[63], 0.4, eps);
  EXPECT_NEAR(sol_vec[64], 0.2, eps);
  EXPECT_NEAR(sol_vec[65], 0, eps);
  EXPECT_NEAR(sol_vec[66], 0, eps);
  EXPECT_NEAR(sol_vec[67], 0, eps);
  EXPECT_NEAR(sol_vec[68], 0.2, eps);
  EXPECT_NEAR(sol_vec[69], 0.4, eps);
  EXPECT_NEAR(sol_vec[70], 0.6, eps);
  EXPECT_NEAR(sol_vec[71], 0.8, eps);
  EXPECT_NEAR(sol_vec[72], 0.8, eps);
  EXPECT_NEAR(sol_vec[73], 0.6, eps);
  EXPECT_NEAR(sol_vec[74], 0.4, eps);
  EXPECT_NEAR(sol_vec[75], 0.2, eps);
  EXPECT_NEAR(sol_vec[76], 0, eps);
  EXPECT_NEAR(sol_vec[77], 0, eps);
  EXPECT_NEAR(sol_vec[78], 0, eps);
  EXPECT_NEAR(sol_vec[79], 0, eps);
  EXPECT_NEAR(sol_vec[80], 0.2, eps);
  EXPECT_NEAR(sol_vec[81], 0.4, eps);
  EXPECT_NEAR(sol_vec[82], 0.6, eps);
  EXPECT_NEAR(sol_vec[83], 0.6, eps);
  EXPECT_NEAR(sol_vec[84], 0.6, eps);
  EXPECT_NEAR(sol_vec[85], 0.4, eps);
  EXPECT_NEAR(sol_vec[86], 0.2, eps);
  EXPECT_NEAR(sol_vec[87], 0, eps);
  EXPECT_NEAR(sol_vec[88], 0, eps);
  EXPECT_NEAR(sol_vec[89], 0, eps);
  EXPECT_NEAR(sol_vec[90], 0, eps);
  EXPECT_NEAR(sol_vec[91], 0, eps);
  EXPECT_NEAR(sol_vec[92], 0.2, eps);
  EXPECT_NEAR(sol_vec[93], 0.4, eps);
  EXPECT_NEAR(sol_vec[94], 0.4, eps);
  EXPECT_NEAR(sol_vec[95], 0.4, eps);
  EXPECT_NEAR(sol_vec[96], 0.4, eps);
  EXPECT_NEAR(sol_vec[97], 0.2, eps);
  EXPECT_NEAR(sol_vec[98], 0, eps);
  EXPECT_NEAR(sol_vec[99], 0, eps);
  EXPECT_NEAR(sol_vec[100], 0, eps);
  EXPECT_NEAR(sol_vec[101], 0, eps);
  EXPECT_NEAR(sol_vec[102], 0, eps);
  EXPECT_NEAR(sol_vec[103], 0, eps);
  EXPECT_NEAR(sol_vec[104], 0.2, eps);
  EXPECT_NEAR(sol_vec[105], 0.2, eps);
  EXPECT_NEAR(sol_vec[106], 0.2, eps);
  EXPECT_NEAR(sol_vec[107], 0.2, eps);
  EXPECT_NEAR(sol_vec[108], 0.2, eps);
  EXPECT_NEAR(sol_vec[109], 0, eps);
  EXPECT_NEAR(sol_vec[110], 0, eps);
  EXPECT_NEAR(sol_vec[111], 0, eps);
  EXPECT_NEAR(sol_vec[112], 0, eps);
  EXPECT_NEAR(sol_vec[113], 0, eps);
  EXPECT_NEAR(sol_vec[114], 0, eps);
  EXPECT_NEAR(sol_vec[115], 0, eps);
  EXPECT_NEAR(sol_vec[116], 0, eps);
  EXPECT_NEAR(sol_vec[117], 0, eps);
  EXPECT_NEAR(sol_vec[118], 0, eps);
  EXPECT_NEAR(sol_vec[119], 0, eps);
  EXPECT_NEAR(sol_vec[120], 0, eps);

  /* END_SOLUTION */
}
