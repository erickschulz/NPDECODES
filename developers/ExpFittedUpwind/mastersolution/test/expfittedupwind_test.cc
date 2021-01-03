#include <gtest/gtest.h>

#include <string>

#include <Eigen/Core>

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>
#include <lf/geometry/geometry.h>

#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>


#include "../expfittedupwind.h"

namespace ExpFittedUpwind::test {

TEST(Bernoulli, Monotinicity) {
    std::vector<double> values = {3.0 , 1.5, 0.5 ,1.0E-5, 1.0E-7, 1.0E-9, -1.0E-9, -1.0E-7, -1.0E-5, -1.0, -5.0, -6.0};

    for(int i = 0; i< values.size()-1; ++i){
        EXPECT_LE(Bernoulli(values[i]), Bernoulli(values[i+1]));
    }
}

TEST(Bernoulli, Cancellation){
    for(int i = 8.0; i >= 10E-12; i/=2.0){
        EXPECT_LE(Bernoulli(i), 2.0);
    }
}

TEST(CompBeta, ConstantPSI){
    auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
    auto mf_Psi = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d /*x*/) {return 3.1;});
	auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    Eigen::VectorXd mu = lf::uscalfe::NodalProjection(*fe_space, mf_Psi);
    double ref = std::exp(3.1);

    auto beta = CompBeta(mesh_p, mu);

    for(auto entity: mesh_p->Entities(1)){
        EXPECT_DOUBLE_EQ((*beta)(*entity), ref);
    }

}

TEST(CompBeta, linearPSI){
     auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
     auto mf_Psi = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) {return 1.0 + x(0) + 2*x(1);});
	 auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
     Eigen::VectorXd mu = lf::uscalfe::NodalProjection(*fe_space, mf_Psi);

    auto beta = CompBeta(mesh_p, mu);
    auto edges = mesh_p->Entities(1);

    //edges of triangle 0 
    EXPECT_DOUBLE_EQ((*beta)(*(edges[0])), std::exp(2.5)*Bernoulli(1.5));
    EXPECT_DOUBLE_EQ((*beta)(*(edges[1])), std::exp(4.0)*Bernoulli(3.0));
    EXPECT_DOUBLE_EQ((*beta)(*(edges[2])), std::exp(4.0)*Bernoulli(3.0));

    //edges of triangle 6:
    EXPECT_DOUBLE_EQ((*beta)(*(edges[12])), std::exp(6.0)*Bernoulli(1.0));
    EXPECT_DOUBLE_EQ((*beta)(*(edges[13])), std::exp(6.0)*Bernoulli(1.0));
    EXPECT_DOUBLE_EQ((*beta)(*(edges[14])), std::exp(6.0)*Bernoulli(0.0));

    //edges of triangle 11
    EXPECT_DOUBLE_EQ((*beta)(*(edges[18])), std::exp(6.0)*Bernoulli(0.0));
    EXPECT_DOUBLE_EQ((*beta)(*(edges[20])), std::exp(6.0)*Bernoulli(-2.5));
    EXPECT_DOUBLE_EQ((*beta)(*(edges[23])), std::exp(8.5)*Bernoulli(2.5));


}


}  // namespace OutputImpedanceBVP::test

