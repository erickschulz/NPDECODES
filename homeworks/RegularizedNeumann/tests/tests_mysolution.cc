#include <gtest/gtest.h>
#include "../mysolution/regNeumann.h"
#include <fstream>
#include <iostream>
#include <lf/io/io.h>

//Test for sub-exercise c with functions h and f being constant
TEST(RegularizedNeumann, solution_test_dropDof_const) {

    //Generate hard-coded mesh file
    std::ofstream myfile;
    myfile.open ("test.msh");
    myfile <<
           "$MeshFormat\n"
           "2.2 0 8\n"
           "$EndMeshFormat\n"
           "$Nodes\n"
           "5\n"
           "1 0 0 0\n"
           "2 1 0 0\n"
           "3 1 1 0\n"
           "4 0 1 0\n"
           "5 0.5 0.5 0\n"
           "$EndNodes\n"
           "$Elements\n"
           "8\n"
           "1 1 2 5 1 1 2\n"
           "2 1 2 5 2 2 3\n"
           "3 1 2 5 3 3 4\n"
           "4 1 2 5 4 4 1\n"
           "5 2 2 8 7 1 2 5\n"
           "6 2 2 8 7 1 5 4\n"
           "7 2 2 8 7 2 3 5\n"
           "8 2 2 8 7 3 4 5\n"
           "$EndElements\n";
    myfile.close();

    //read in hard-coded mesh file
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader reader(std::move(mesh_factory), "test.msh");
    auto mesh_p = reader.mesh();

    //Delete the hard-coded mesh file
    std::remove("test.msh");

    //source and boundary functions for testing
    const auto f = lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> double { return 1.0; });
    const auto h = lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> double { return 1.0; });

    auto fe_space =
            std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

    // Compute solution
    auto result_c = RegularizedNeumann::getGalerkinLSE_dropDof(fe_space, f, h);


    //Now we compare the returned values to hard coded results
    const double eps = 1e-10;

    //Check Matrix
    Eigen::MatrixXd solution_mat_c (5,5);
    solution_mat_c << 1, 0, 0, 0, 0,
            0, 1, 0, 0, -1,
            0, 0, 1, 0, -1,
            0, 0, 0, 1, -1,
            0, -1, -1,-1,4;
    //Compare with expected results
    for (int i = 0; i <5 ; ++i) {
        for (int j = 0; j < 5; ++j) {
            EXPECT_NEAR(solution_mat_c(i, j), result_c.first.coeff(i, j), eps);
        }
    }
    //Check rhs vector
    Eigen::VectorXd solution_vec_c (5);
    solution_vec_c << 0 , 1.1666666667, 1.1666666667, 1.1666666667, 0.33333333333;
    //Compare with expected results
    for (int i = 0; i <solution_vec_c.size() ; ++i) {
        EXPECT_NEAR(solution_vec_c(i),result_c.second(i) , eps);
    }
}

//Test for sub-exercise c with functions h and f not being constant
TEST(RegularizedNeumann, solution_test_dropDof_gen) {

    //Generate hard-coded mesh file
    std::ofstream myfile;
    myfile.open ("test.msh");
    myfile <<
           "$MeshFormat\n"
           "2.2 0 8\n"
           "$EndMeshFormat\n"
           "$Nodes\n"
           "5\n"
           "1 0 0 0\n"
           "2 1 0 0\n"
           "3 1 1 0\n"
           "4 0 1 0\n"
           "5 0.5 0.5 0\n"
           "$EndNodes\n"
           "$Elements\n"
           "8\n"
           "1 1 2 5 1 1 2\n"
           "2 1 2 5 2 2 3\n"
           "3 1 2 5 3 3 4\n"
           "4 1 2 5 4 4 1\n"
           "5 2 2 8 7 1 2 5\n"
           "6 2 2 8 7 1 5 4\n"
           "7 2 2 8 7 2 3 5\n"
           "8 2 2 8 7 3 4 5\n"
           "$EndElements\n";
    myfile.close();

    //read in hard-coded mesh file
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader reader(std::move(mesh_factory), "test.msh");
    auto mesh_p = reader.mesh();

    //Delete the hard-coded mesh file
    std::remove("test.msh");

    //source and boundary functions for testing
    const auto f = lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> double { return x(0)+x(1); });
    const auto h = lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> double { return x(0)+x(1); });

    auto fe_space =
            std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

    // Compute solution
    auto result_c = RegularizedNeumann::getGalerkinLSE_dropDof(fe_space, f, h);

    //Now we compare the returned values to hard coded results
    const double eps = 1e-5;

    //Check Matrix
    Eigen::MatrixXd solution_mat_c (5,5);
    solution_mat_c << 1, 0, 0, 0, 0,
            0, 1, 0, 0, -1,
            0, 0, 1, 0, -1,
            0, 0, 0, 1, -1,
            0, -1, -1,-1,4;
    //Compare with expected results
    for (int i = 0; i <5 ; ++i) {
        for (int j = 0; j < 5; ++j) {
            EXPECT_NEAR(solution_mat_c(i, j), result_c.first.coeff(i, j), eps);
        }
    }
    //Check rhs vector
    Eigen::VectorXd solution_vec_c (5);
    solution_vec_c << 0 , 1.1666666667, 1.91667, 1.1666666667, 0.33333333333;
    //Compare with expected results
    for (int i = 0; i <solution_vec_c.size() ; ++i) {
        EXPECT_NEAR(solution_vec_c(i),result_c.second(i) , eps);
    }
}

//Test for sub-exercise f with functions h and f being constant
TEST(RegularizedNeumann, solution_test_augment_const) {

    //Generate hard-coded mesh file
    std::ofstream myfile;
    myfile.open ("test.msh");
    myfile <<
           "$MeshFormat\n"
           "2.2 0 8\n"
           "$EndMeshFormat\n"
           "$Nodes\n"
           "5\n"
           "1 0 0 0\n"
           "2 1 0 0\n"
           "3 1 1 0\n"
           "4 0 1 0\n"
           "5 0.5 0.5 0\n"
           "$EndNodes\n"
           "$Elements\n"
           "8\n"
           "1 1 2 5 1 1 2\n"
           "2 1 2 5 2 2 3\n"
           "3 1 2 5 3 3 4\n"
           "4 1 2 5 4 4 1\n"
           "5 2 2 8 7 1 2 5\n"
           "6 2 2 8 7 1 5 4\n"
           "7 2 2 8 7 2 3 5\n"
           "8 2 2 8 7 3 4 5\n"
           "$EndElements\n";
    myfile.close();

    //read in hard-coded mesh file
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader reader(std::move(mesh_factory), "test.msh");
    auto mesh_p = reader.mesh();

    //Delete the hard-coded mesh file
    std::remove("test.msh");

    //source and boundary functions for testing
    const auto f = lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> double { return 1.0; });
    const auto h = lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> double { return 1.0; });

    auto fe_space =
            std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

    // Compute solution
    auto result_f = RegularizedNeumann::getGalerkinLSE_augment(fe_space, f, h);

    //Compare the returned values to hard coded results
    const double eps = 1e-5;

    //Check Matrix
    Eigen::MatrixXd solution_mat_f (6,6);
    solution_mat_f << 1, 0, 0, 0, -1, 0.1666667,
            0, 1, 0, 0, -1, 0.1666667,
            0, 0, 1, 0, -1, 0.1666667,
            0, 0, 0, 1, -1, 0.1666667,
            -1, -1, -1,-1,4, 0.3333333,
            0.166667, 0.166667, 0.166667, 0.166667, 0.3333333, 0;
    //Compare with expected results
    for (int i = 0; i <6 ; ++i) {
        for (int j = 0; j < 6; ++j) {
            EXPECT_NEAR(solution_mat_f(i, j), result_f.first.coeff(i, j), eps);
        }
    }

    //Check rhs vector
    Eigen::VectorXd solution_vec_f (6);
    solution_vec_f <<  1.1666666667, 1.1666666667, 1.1666666667, 1.1666666667, 0.33333333333, 0;

    for (int i = 0; i <solution_vec_f.size() ; ++i) {
        EXPECT_NEAR(solution_vec_f(i),result_f.second(i) , eps);
    }

}

//Test for sub-exercise f with functions h and f not being constant
TEST(RegularizedNeumann, solution_test_augment_gen) {

    //Generate hard-coded mesh file
    std::ofstream myfile;
    myfile.open ("test.msh");
    myfile <<
           "$MeshFormat\n"
           "2.2 0 8\n"
           "$EndMeshFormat\n"
           "$Nodes\n"
           "5\n"
           "1 0 0 0\n"
           "2 1 0 0\n"
           "3 1 1 0\n"
           "4 0 1 0\n"
           "5 0.5 0.5 0\n"
           "$EndNodes\n"
           "$Elements\n"
           "8\n"
           "1 1 2 5 1 1 2\n"
           "2 1 2 5 2 2 3\n"
           "3 1 2 5 3 3 4\n"
           "4 1 2 5 4 4 1\n"
           "5 2 2 8 7 1 2 5\n"
           "6 2 2 8 7 1 5 4\n"
           "7 2 2 8 7 2 3 5\n"
           "8 2 2 8 7 3 4 5\n"
           "$EndElements\n";
    myfile.close();

    //read in hard-coded mesh file
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader reader(std::move(mesh_factory), "test.msh");
    auto mesh_p = reader.mesh();

    //Delete the hard-coded mesh file
    std::remove("test.msh");

    //source and boundary functions for testing
    const auto f = lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> double { return x(0)+x(1); });
    const auto h = lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> double { return x(0)+x(1); });

    auto fe_space =
            std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

    // Compute solution
    auto result_f = RegularizedNeumann::getGalerkinLSE_augment(fe_space, f, h);

    //Compare the returned values to hard coded results
    const double eps = 1e-5;

    //Check Matrix
    Eigen::MatrixXd solution_mat_f (6,6);
    solution_mat_f << 1, 0, 0, 0, -1, 0.1666667,
            0, 1, 0, 0, -1, 0.1666667,
            0, 0, 1, 0, -1, 0.1666667,
            0, 0, 0, 1, -1, 0.1666667,
            -1, -1, -1,-1,4, 0.3333333,
            0.166667, 0.166667, 0.166667, 0.166667, 0.3333333, 0;
    //Compare with expected results
    for (int i = 0; i <6 ; ++i) {
        for (int j = 0; j < 6; ++j) {
            EXPECT_NEAR(solution_mat_f(i, j), result_f.first.coeff(i, j), eps);
        }
    }

    //Check rhs vector
    Eigen::VectorXd solution_vec_f (6);
    solution_vec_f <<  0.416667, 1.1666666667, 1.91667, 1.1666666667, 0.33333333333, 0;

    for (int i = 0; i <solution_vec_f.size() ; ++i) {
        EXPECT_NEAR(solution_vec_f(i),result_f.second(i) , eps);
    }
}

