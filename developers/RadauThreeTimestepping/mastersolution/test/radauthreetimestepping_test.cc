/**
 * @file radauthreetimestepping_test.cc
 * @brief NPDE homework "RadauThreeTimestepping" code
 * @author Tobias Rohner
 * @date 16.03.2020
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>

#include <Eigen/Core>

#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/geometry/geometry.h>

#include "../radauthreetimestepping.h"

namespace RadauThreeTimestepping::test {

TEST(RadauThreeTimestepping, TrapRuleLinFEElemVecProvider) {
    // Get some triangular test mesh
    auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
    // Define some easy functions to test the provider with
    auto f1 = [](const Eigen::Vector2d &x) -> double { return 0; };
    auto f2 = [](const Eigen::Vector2d &x) -> double { return 1; };
    auto f3 = [](const Eigen::Vector2d &x) -> double { return x[0]; };
    // Check the element vector for each triangle
    RadauThreeTimestepping::TrapRuleLinFEElemVecProvider f1p(f1);
    RadauThreeTimestepping::TrapRuleLinFEElemVecProvider f2p(f2);
    RadauThreeTimestepping::TrapRuleLinFEElemVecProvider f3p(f3);
    for (const auto tria : mesh_p->Entities(0)) {
	const auto geom = tria->Geometry();
	auto ev1 = f1p.Eval(*tria);
	auto ev2 = f2p.Eval(*tria);
	auto ev3 = f3p.Eval(*tria);
	ASSERT_TRUE(ev1.isApprox(Eigen::Vector3d::Zero()));
	ASSERT_TRUE(ev2.isApprox(Eigen::Vector3d::Constant(lf::geometry::Volume(*geom)/3)));
	Eigen::Vector3d e3_correct = lf::geometry::Volume(*geom)/3 * lf::geometry::Corners(*geom).row(0);
	ASSERT_TRUE(ev3.isApprox(e3_correct));
    }
}

}   // end namespace RadauThreeTimestepping::test
