/**
 * @file advectionfv2d.cc
 * @brief NPDE homework AdvectionFV2D code
 * @author Philipp Egg
 * @date 21.06.2020
 * @copyright Developed at ETH Zurich
 */

#include "advectionfv2d.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

namespace AdvectionFV2D {

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd dummyFunction(double x, int n) {
#if SOLUTION
  // Appears only in mastersolution
  return Eigen::Vector2d::Constant(1.0);
#else
  // Appears only in mysolution and templates
  return Eigen::Vector2d::Zero();
#endif
}

// A helper function to compute the normal vectors
Eigen::Matrix<double, 2, 3> gradbarycoordinates(
    const Eigen::Matrix<double, 2, 3> &triangle){
  Eigen::Matrix3d X;
  
  // solve for the coefficients of the barycentric coordinate functions
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = triangle.transpose();
  return X.inverse().block<2, 3>(1, 0);
}

std::shared_ptr<lf::mesh::utils::CodimMeshDataSet
    <Eigen::Matrix<double, 2, Eigen::Dynamic>>> 
    computeCellNormals(std::shared_ptr<const lf::mesh::Mesh> mesh_p){

#if SOLUTION
  // Appears only in mastersolution

  // Initialize datastruture for the result
  lf::mesh::utils::CodimMeshDataSet
      <Eigen::Matrix<double,2,Eigen::Dynamic>> result(mesh_p, 0);

  // Compute normal vectors
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
    const lf::geometry::Geometry *geo_p = cell->Geometry();
    const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

    if (corners.cols() == 3){
      Eigen::Matrix<double, 2, Eigen::Dynamic> normal_vectors(2, 3);

      // Compute normal vectors
      Eigen::Matrix<double, 2, 3> bary_coord = gradbarycoordinates(corners);

      // Reorder computed normal vectors
      // normal_vectors[0] -> normal_vector of edge[0]
      normal_vectors.col(0) = -(bary_coord.col(2)).normalized();
      normal_vectors.col(1) = -(bary_coord.col(0)).normalized();
      normal_vectors.col(2) = -(bary_coord.col(1)).normalized();

      // Save the matrix in our datastructure
      result(*cell) = normal_vectors;
    } else { // corners.cols() == 4
      Eigen::Matrix<double, 2, Eigen::Dynamic> normal_vectors(2, 4);

      // Split the quadrilateral into two triangles (1,2,3) and (1,3,4)
      Eigen::Matrix<double, 2, 3> tria_1; 
      Eigen::Matrix<double, 2, 3> tria_2; 
      tria_1 << corners.col(0), corners.col(1), corners.col(2);
      tria_2 << corners.col(0), corners.col(2), corners.col(3);

      // Compute normal vectors
      Eigen::Matrix<double, 2, 3> bary_coord_1 = gradbarycoordinates(tria_1);
      Eigen::Matrix<double, 2, 3> bary_coord_2 = gradbarycoordinates(tria_2);

      // Reorder computed normal vectors
      // normal_vectors[0] -> normal_vector of edge[0]
      normal_vectors.col(0) = -(bary_coord_1.col(2)).normalized();
      normal_vectors.col(1) = -(bary_coord_1.col(0)).normalized();
      normal_vectors.col(2) = -(bary_coord_2.col(0)).normalized();
      normal_vectors.col(3) = -(bary_coord_2.col(1)).normalized();

      // Save the matrix in our datastructure
      result(*cell) = normal_vectors;
    }
  }
  return std::make_shared<lf::mesh::utils::CodimMeshDataSet
      <Eigen::Matrix<double,2,Eigen::Dynamic>>>(result);
  
#else
  // Appears only in mysolution and templates
#endif
}

std::shared_ptr<lf::mesh::utils::CodimMeshDataSet
    <std::array<const lf::mesh::Entity*, 4>>> getAdjacentCellPointers(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p){
  // Initialize auxilary object
  lf::mesh::utils::CodimMeshDataSet
      <std::array<const lf::mesh::Entity*, 2>> 
      aux_obj(mesh_p, 1, {nullptr, nullptr});

  // Iterate over every cell
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)){
    auto cell_edges = cell->SubEntities(1);

    // Iterate over every edge of the cell
    for (const lf::mesh::Entity *edge : cell_edges){
      // If aux_obj at index 0 was not set, save cell there
      // otherwise save at the second position
      // (The first position has to be set; the second might be set)
      if (aux_obj(*edge)[0] == nullptr){
        aux_obj(*edge)[0] = cell;
      } else if (aux_obj(*edge)[1] == nullptr){
        aux_obj(*edge)[1] = cell;
      } else {
        throw std::runtime_error("Error in aux_obj");
      }
    }
  }

  // Initialize datastructure for result
  lf::mesh::utils::CodimMeshDataSet
      <std::array<const lf::mesh::Entity*, 4>> result(
      mesh_p, 0, {nullptr, nullptr, nullptr, nullptr});

  // Collect the objects of the auxilary object
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)){
    auto cell_edges = cell->SubEntities(1);
    int counter = 0;
    for (const lf::mesh::Entity *edge : cell_edges){
      if (aux_obj(*edge)[0] != cell){
        result(*cell)[counter] = aux_obj(*edge)[0];
      } else{
        result(*cell)[counter] = aux_obj(*edge)[1];
      }
      counter++;
    }
  }

  return std::make_shared <lf::mesh::utils::CodimMeshDataSet
      <std::array<const lf::mesh::Entity*, 4>>>(result);
}

// Function returning the barycenter of TRIA or QUAD
Eigen::Vector2d barycenter(const Eigen::MatrixXd corners){
  Eigen::Vector2d midpoint;
  if (corners.cols() == 3){
    midpoint = (corners.col(0) + corners.col(1) + corners.col(2)) / 3.0;
  } else if (corners.cols() == 4) {
    midpoint = (corners.col(0) + corners.col(1) + corners.col(2) + 
                corners.col(3)) / 4.0;
  } else {
    throw std::runtime_error("Wrong geometrie in barycenter()");
  }
  return midpoint;
}

double computeHmin(std::shared_ptr<const lf::mesh::Mesh> mesh_p){
  //Vector to store the distances between cells
  std::vector<double> min_h;
  
  // Get Adjectent Cells
  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet
      <std::array<const lf::mesh::Entity*, 4>>>
      adjacentCells = AdvectionFV2D::getAdjacentCellPointers(mesh_p);

  //Iterate over all cells
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
    const lf::geometry::Geometry* geo_p = cell->Geometry();
    const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

    //Compute the barycenter of the cell
    Eigen::Vector2d cur_midpoint = barycenter(corners);

    //Iterate over all adjecent cells
    for (const lf::mesh::Entity *neighbour_cell: (*adjacentCells)(*cell)){
      //Check that the neighbor exists
      if (neighbour_cell != nullptr){
        const lf::geometry::Geometry *geo_p_neighbour = 
            neighbour_cell->Geometry();
        const Eigen::MatrixXd neighbour_corners = 
            lf::geometry::Corners(*geo_p_neighbour);

        //Compute barycenter of neighbour cell
        Eigen::Vector2d neighbour_midpoint = 
            barycenter(neighbour_corners);

        //Compute distances between cells
        Eigen::Vector2d diff = cur_midpoint - neighbour_midpoint;
        double distance = diff.norm();

        //Store value in vector
        min_h.push_back(distance);
      }
    }
  }

  //Find the minimum value in the vector and return it
  double result = *std::min_element(std::begin(min_h), std::end(min_h));
  return result;
}

// TODO: Problem description Vector2d?? -> changed to VectorXd
Eigen::VectorXd simulateAdvection(const lf::assemble::DofHandler &dofh){
  //Specifiy final time
  double T_max = 1.0;

  Eigen::Vector2d x0(0.8, 0.2);
  //////////////////////////////////////////////////////////////////////
  //TODO Inconsistency (x - x0) -> ||x - x0||
  //////////////////////////////////////////////////////////////////////
  // Set up lambda function u0 and beta
  double d = 0.2;
  auto u0 = [x0, d](Eigen::Vector2d x) -> double {
    double dist = (x - x0).norm();
    if (dist < d){
      return std::pow(std::cos(M_PI / (2 * d) * dist), 2);
    } else {
      return 0.0;
    }
  };
  auto beta = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x[1], x[0]) / std::sqrt(2.0);
  };

  // Get number of dofs and
  // inititialize a vector for the initial condition
  int num_dof = dofh.NumDofs();
  Eigen::VectorXd u0_h(num_dof);

  // Iterate over all cells
  for (const lf::mesh::Entity *cell : dofh.Mesh()->Entities(0)) {
    const lf::geometry::Geometry *geo_p = cell->Geometry();
    const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

    // Compute barycenter of the cell
    Eigen::Vector2d x = barycenter(corners);

    // Find the index of the DOF for the cell
    int idx = dofh.GlobalDofIndices(*cell)[0];

    double distance = (x - x0).norm();
    if (distance < d){
      u0_h[idx] = u0(x);
    } else {
      u0_h[idx] = 0.0;
    }
  }

  // Set up the number of steps according to the CFL condition
  int M = int ((T_max / computeHmin(dofh.Mesh())) + 2);

  // Run simulation
  Eigen::VectorXd result = solveAdvection2D(dofh, beta, u0_h, T_max , M);

  return result;
}

Eigen::VectorXd refSolution(const lf::assemble::DofHandler &dofh){
  // Specifiy final time
  double T_max = 1.0;

  Eigen::Vector2d x0(0.8, 0.2);
  // Set up lambda function u0 and beta
  double d = 0.2;
  auto u0 = [x0, d](Eigen::Vector2d x) -> double {
    double dist = (x - x0).norm();
    if (dist < d){
      return std::pow(std::cos(M_PI / (2 * d) * dist), 2);
    } else {
      return 0.0;
    }
  };

  // Setup inverted phi^-1  
  Eigen::Matrix2d phi_inv;
  double t_mod = T_max / std::sqrt(2.0);
  phi_inv << std::cos(t_mod), std::sin(t_mod), 
            -std::sin(t_mod), std::cos(t_mod);

  // Initialize vector for the result
  int num_dof = dofh.NumDofs();
  Eigen::VectorXd ref_solution(num_dof);

  // Iterate over all cells
  for (const lf::mesh::Entity *cell : dofh.Mesh()->Entities(0)) {
    const lf::geometry::Geometry *geo_p = cell->Geometry();
    const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

    // Compute the barycenter of the cell
    Eigen::Vector2d x = barycenter(corners);

    // Find the index of the DOF for the cell
    int idx = dofh.GlobalDofIndices(*cell)[0];

    // Compute exact solution
    Eigen::Vector2d phi_inv_x = phi_inv * x;
    ref_solution[idx] =  u0(phi_inv_x);
  }  
  return ref_solution;
}
/* SAM_LISTING_END_1 */

}  // namespace AdvectionFV2D
