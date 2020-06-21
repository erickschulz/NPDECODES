/**
 * @file advectionfv2d_main.cc
 * @brief NPDE homework AdvectionFV2D code
 * @author Philipp Egg
 * @date 21.06.2020
 * @copyright Developed at ETH Zurich
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>
#include <lf/refinement/refinement.h>

#include "advectionfv2d.h"

int main() {
	
  //TODO inconsistancy: g vs. G
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1./3.);
  
  auto beta = [](Eigen::Vector2d x) -> Eigen::Vector2d { 
      return Eigen::Vector2d(-x[1], x[0]); };
  
  // Perferforming some manual tests for computeCellNormals and 
  // getAdjacentCellPointers (Tests were OK)
  /*
  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet
      <Eigen::Matrix<double,2,Eigen::Dynamic>>> 
      normal_vectors = AdvectionFV2D::computeCellNormals(mesh_p);

  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet
      <std::array<const lf::mesh::Entity*,4>>>
      adjacentCells = AdvectionFV2D::getAdjacentCellPointers(mesh_p);

  int element_idx = 0;
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
    const lf::geometry::Geometry* geo_p = cell->Geometry();
    const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);
    
    Eigen::Vector2d cur_midpoint = AdvectionFV2D::barycenter(corners);
    
    std::cout << "Element NR: " << element_idx << std::endl;
    std::cout << "Corners: \n" << corners << std::endl;
    std::cout << "Barycenter: \n" << cur_midpoint << std::endl;
    std::cout << "Normal Vectors\n" << (*normal_vectors)(*cell) << std::endl;
    
    for (const lf::mesh::Entity* neighbour_cell: (*adjacentCells)(*cell)){
		if (neighbour_cell != nullptr){
		  int neighbour_ind = 0;
          for (const lf::mesh::Entity *cell2 : mesh_p->Entities(0)){
			  if (neighbour_cell == cell2){
				  std::cout << "Neighbour at: " << neighbour_ind << std::endl;
			  }
			  neighbour_ind++;
		  }
		} else {
			std::cout << "Neighbour at: - "<< std::endl;
		}
	}
    element_idx++;
  }
  */
  
  
  // Output B_Matrix
  /*
  const lf::assemble::UniformFEDofHandler dofh(mesh_p, 
      {{lf::base::RefEl::kPoint(), 0},
      {lf::base::RefEl::kSegment(), 0},
      {lf::base::RefEl::kTria(), 1},
      {lf::base::RefEl::kQuad(), 1}});
  
  Eigen::SparseMatrix<double> B_Matrix = 
      AdvectionFV2D::initializeMOLODEMatrix(dofh, beta);
      
  std::cout << "B_Matrix\n" << B_Matrix << std::endl;
  */

  //Testing Output minH
  /*
  double minH = AdvectionFV2D::computeHmin(mesh_p);
  std::cout << "minh " << minH << std::endl;
  */
  
  //Test of Task 8-8.n
  //Eigen::VectorXd result = AdvectionFV2D::simulateAdvection(dofh);

  // Task 8-8.o
  //////////////////////////////////////////////////////////////////////
  // TODO: Change to refinement_level = 5
  //////////////////////////////////////////////////////////////////////
  // Generate a mesh hierarchy
  auto mesh_seq_p{ 
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p,4)};
  int num_meshes = mesh_seq_p->NumLevels();
  
  int array_num_cells[num_meshes] = {};
  double array_l2error[num_meshes] = {};
  
  // Iterate over all mesh levels
  for (int level = 0; level < num_meshes; level++){
    std::cout << "Computing L2Error for level: " << level << std::endl;
    
    // Get the current mesh
    auto cur_mesh = mesh_seq_p->getMesh(level);
	  
    // Create a DOF Hander for the current mesh
    const lf::assemble::UniformFEDofHandler cur_dofh(cur_mesh, 
        {{lf::base::RefEl::kPoint(), 0},
        {lf::base::RefEl::kSegment(), 0},
        {lf::base::RefEl::kTria(), 1},
        {lf::base::RefEl::kQuad(), 1}});
	
    //Get Result from Simulation
    Eigen::VectorXd result = AdvectionFV2D::simulateAdvection(cur_dofh);
    
    //Get Exact Result
    Eigen::VectorXd ref_solution = AdvectionFV2D::refSolution(cur_dofh);
    
    //Compute L2 Error
    double l2_error = 0;
    for (const lf::mesh::Entity *cell : cur_mesh->Entities(0)) {
      const lf::geometry::Geometry* geo_p = cell->Geometry();
      double area = lf::geometry::Volume(*geo_p);
      int idx = cur_dofh.InteriorGlobalDofIndices(*cell)[0];
      l2_error += std::pow( (result[idx] - ref_solution[idx]), 2) * area;
    }
    
    // Store the results in the arrays
    array_num_cells[level] = cur_dofh.NumDofs();
    array_l2error[level] = l2_error;
  }


  //Task 8-8.q
  double array_clf_thres[num_meshes] = {};
  
  //Iterate over all mesh levels
  for (int level = 0; level < num_meshes; level++){
    // Get the current mesh
    auto cur_mesh = mesh_seq_p->getMesh(level);
    
    // Create a DOF Hander for the current mesh
    const lf::assemble::UniformFEDofHandler cur_dofh(cur_mesh, 
        {{lf::base::RefEl::kPoint(), 0},
        {lf::base::RefEl::kSegment(), 0},
        {lf::base::RefEl::kTria(), 1},
        {lf::base::RefEl::kQuad(), 1}});
    
    //Save the CFL threshold to array for T = 1.0
    std::cout << "Computing threshold for level: " << level << std::endl;
    array_clf_thres[level] = 
        AdvectionFV2D::findCFLthreshold(cur_dofh, beta, 1.0);
  }
  
  //Write Output file of Task 8-8.o and 8-8.q
  std::ofstream csv_file;
  csv_file.open ("advectionfv2d.csv");
  for (int level = 0; level < num_meshes; level++){
    std::cout << "Cells: "<< array_num_cells[level] 
              << " | Error: " << array_l2error[level] 
              << std::endl;
    csv_file << array_num_cells[level] << "," 
             << array_l2error[level] << "," 
             << array_clf_thres[level] << "\n";
  }
  csv_file.close();

  //const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
  //    Eigen::DontAlignCols, ", ", "\n");
  //std::cout << v.transpose().format(CSVFormat) << std::endl;
  return 0;
}
