#include <iostream>
// Eigen includes
#include <Eigen/Core>
#include <Eigen/SparseLU>
// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

#include "../waveabc2d.h"

#include <gtest/gtest.h>

namespace WaveABC2D::test {
TEST(WaveABC2D, scalarImplicitTimestepping) {
    double eps = 1.0e-5;
    double epsilon = 0.5;
    int M = 10;
    Eigen::VectorXd test_sol_vec(M+1);
    test_sol_vec <<1.0,
                1.08531,
                1.14196,
                1.17159,
                1.1762,
                1.15807,
                1.11971,
                1.06376,
                0.992949,
                0.910024,
                0.817708;

    Eigen::VectorXd tested_solution = scalarImplicitTimestepping(epsilon,M);

    for (int i = 0; i<M+1;i++){
        ASSERT_NEAR(tested_solution(i), test_sol_vec(i),eps);
    }
}

//Test with constant initial data
TEST(WaveABC2D,WaveABC2DTimestepper_const){
    auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
    auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    auto rho = [](Eigen::Vector2d) -> double { return 4.0; };
    auto mu0 = [](const Eigen::Vector2d &x) -> double {
        return 1.0;
    };
    auto nu0 = [](const Eigen::Vector2d &x) -> double { 
        return 1.0; 
    };
    auto stepper = WaveABC2DTimestepper<decltype(rho), decltype(mu0), decltype(nu0)>(fe_space_p, rho, 500, 1.0);
    auto tested_solution =  stepper.solveWaveABC2D(mu0, nu0);
    Eigen::VectorXd reference_solution(13);
    reference_solution <<  1.57839,
                    1.71067,
                    1.57317,
                    2.11469,
                    2.14686,
                    1.68006,
                    1.71359,
                    2.15184,
                    2.11977,
                    1.64308,
                    1.57195,
                    1.74435,
                    1.62919;
    ASSERT_TRUE(tested_solution.size()==reference_solution.size());
    for(int i = 0; i<tested_solution.size(); i++){
        ASSERT_NEAR(tested_solution(i), reference_solution(i), 1.0e-5);
    }
}

//Test with original initial data
TEST(WaveABC2D,WaveABC2DTimestepper){
    auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
    auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    auto rho = [](Eigen::Vector2d) -> double { return 4.0; };

    auto mu0 = [](const Eigen::Vector2d &x) -> double {
        return std::sin(x.norm());
    };
    auto nu0 = [](const Eigen::Vector2d &x) -> double { 
        return std::cos(x(1)); 
    };
    auto stepper = WaveABC2DTimestepper<decltype(rho), decltype(mu0), decltype(nu0)>(fe_space_p, rho, 500, 1.0);
    auto tested_solution =  stepper.solveWaveABC2D(mu0, nu0);
    Eigen::VectorXd reference_solution(13);
    reference_solution <<  0.97508,
                1.50783,
                0.876215,
                1.47616,
                1.1743,
                0.614315,
                0.882905,
                0.342194,
                0.0585756,
                -0.343122,
                -0.232691,
                -0.622208,
                -1.36191;
    ASSERT_TRUE(tested_solution.size()==reference_solution.size());
    for(int i = 0; i<tested_solution.size(); i++){
        ASSERT_NEAR(tested_solution(i), reference_solution(i), 1.0e-5);
    }
}

}
