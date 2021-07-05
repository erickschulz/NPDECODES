/**
 * @file upwindquadrature_main.cc
 * @brief NPDE homework template main
 * @author Philippe Peter
 * @date June 2020
 * @copyright Developed at SAM, ETH Zurich
 */
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <memory>

#include "standard_galerkin.h"
#include "convection_emp.h"
#include "sample_meshfunction.h"
#include "solve_upwind.h"
#include "streamline_upwind.h"

int main(){
    //parameter functions:
    //boundary conditions
    const auto g = [](const Eigen::Vector2d &x) {
       return x(0) > x(1)? 1.0:0.0;
    };
    // velocity field
    const auto v = [](const Eigen::Vector2d &x) {
        return Eigen::Vector2d(1.0,1.0);
    };
    const auto eps = [](const Eigen::Vector2d &x) {
       return std::pow(10,-10);
    };

    const auto f = [](const Eigen::Vector2d &x) {
       return 0.0;
    };

    //construct mesh:
    int M = 49*2;

     // MESH CONSTRUCTION
    // construct a triangular tensor product mesh on the unit square
/*  
    std::unique_ptr<lf::mesh::MeshFactory> mesh_factory_ptr =
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::mesh::utils::TPTriagMeshBuilder builder(std::move(mesh_factory_ptr));
    builder.setBottomLeftCorner(Eigen::Vector2d{0.0, 0.0})
        .setTopRightCorner(Eigen::Vector2d{1.0, 1.0})
        .setNumXCells(M)
        .setNumYCells(M);
    std::shared_ptr<lf::mesh::Mesh> mesh_p = builder.Build();
*/
    std::string mesh_file = "../../../lecturecodes/ConvectionDiffusion/mesh_square.gmsh";
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader reader(std::move(mesh_factory),mesh_file);
    auto mesh_p = reader.mesh();

    // DOF HANDLER & FINITE ELEMENT SPACE
    // Construct dofhanlder for linear finite elements on the mesh.
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

    Eigen::VectorXd sol_standard = ConvectionDiffusion::SolveCDBVPStandardGalerkin(fe_space,eps,v,f,g);
    lf::fe::MeshFunctionFE sol_standard_mf(fe_space, sol_standard);

    Eigen::VectorXd sol_stable = ConvectionDiffusion::SolveCDBVPUpwind(fe_space,eps,v,f,g);
    lf::fe::MeshFunctionFE sol_upwind_mf(fe_space, sol_stable);

    Eigen::VectorXd sol_supg = ConvectionDiffusion::SolveCDBVPSupg(fe_space,eps,v,f,g);
    lf::fe::MeshFunctionFE sol_supg_mf(fe_space, sol_supg);


    auto gamma = [](double t){ return Eigen::Vector2d(t,1-t);};
    ConvectionDiffusion::SampleMeshFunction("standard_galerkin.txt", mesh_p, gamma, sol_standard_mf, 100);
    ConvectionDiffusion::SampleMeshFunction("upwind.txt", mesh_p, gamma, sol_upwind_mf, 100);
    ConvectionDiffusion::SampleMeshFunction("supg.txt", mesh_p, gamma, sol_supg_mf, 100);

    
    return 0;
}