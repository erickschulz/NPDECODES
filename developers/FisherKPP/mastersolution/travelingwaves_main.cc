/** @file travelingwaves_main.cc
 *  @brief Homework Problem FisherKPP
 *  @author Amélie Justine Loher
 *  @date 12.05.20
 *  @copyright Developed at SAM, ETH Zurich
 */

#include "fisherkpp.cc"

#include <cstdlib>
#include <iostream>
#include <string>
#include <stdexcept>

#include <lf/io/io.h>

using namespace FisherKPP;

int main(int /*argc*/, char ** /*argv*/){

  // Obtain mesh 
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), CURRENT_SOURCE_DIR "/../meshes/island.msh");
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
  // Finite Element Space 
  auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Dofhandler 
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Initial Population density 
  Eigen::VectorXd u0(N_dofs); u0.setZero();
  u0(321) = 0.3;            
  u0(567) = 0.3;

  // Diffusion Coefficient
  auto c = [] (Eigen::Vector2d x) -> double { return 1.2;};

  // Growth Factor
  double lambda = 2.1;

  // Carrying capacity
  Eigen::VectorXd K(N_dofs); K.setZero();
  K = 0.8 * Eigen::VectorXd::Ones(N_dofs);

  // Time Steps
  unsigned int m = 100;
  double T = 1.;

  // Now we may compute the solution.
  StrangSplit StrangSplitter(fe_space, T, m, lambda, c);
  
  Eigen::MatrixXd sol(N_dofs, 5);
  
  std::cout << "You are running the model problem, i.e. you solve the Fisher equation on the model domain. "
  			<< std::endl;

  sol.col(0) = StrangSplitter.Evolution(K, u0);
  std::cout << "Solution after 100 timesteps." << std::endl;
  sol.col(1) = StrangSplitter.Evolution(K, sol.col(0));
  std::cout << "Solution after 200 timesteps." << std::endl;
  sol.col(2) = StrangSplitter.Evolution(K, sol.col(1));
  std::cout << "Solution after 300 timesteps." << std::endl;
  sol.col(3) = StrangSplitter.Evolution(K, sol.col(2));
  std::cout << "Solution after 400 timesteps." << std::endl;
  sol.col(4) = StrangSplitter.Evolution(K, sol.col(3));
  std::cout << "Solution after 500 timesteps." << std::endl;

  // Use VTK-Writer for Visualization of solution.
  #if SOLUTION 
  
  for(int k = 1; k < 6; k++) {
    
	std::stringstream filename;
    filename << "sol" << k << ".vtk";
	
	lf::io::VtkWriter vtk_writer(mesh_p, filename.str());
	auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
    for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
      nodal_data->operator()(dofh.Entity(global_idx)) = sol.col(k-1)[global_idx];
    }

    vtk_writer.WritePointData("sol", *nodal_data);
  }

  #else
  //====================
  // Your code goes here
  //====================
  #endif

  return 0;
 
}
