/**
 * @file advectionfv2d_main.cc
 * @brief NPDE homework AdvectionFV2D code
 * @author Philipp Egg
 * @date 21.06.2020
 * @copyright Developed at ETH Zurich
 */

#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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
#include <lf/quad/quad_rule.h>
#include <lf/refinement/refinement.h>

#include "advectionfv2d.h"

int main() {
  // TODO inconsistancy: g vs. G
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1. / 3.);
  // Or use a tria-mesh
  // auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1. / 3.);

  // Define velocity field
  // Note that the problem description requires ||B|| <= 1
  auto beta = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x[1], x[0]) / std::sqrt(2.0);
  };

  // Perferforming some manual tests for computeCellNormals and
  /*
  // getAdjacentCellPointers (Tests were OK)
  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet
      <Eigen::Matrix<double, 2, Eigen::Dynamic>>>
      normal_vectors = AdvectionFV2D::computeCellNormals(mesh_p);

  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet
      <std::array<const lf::mesh::Entity*, 4>>>
      adjacentCells = AdvectionFV2D::getAdjacentCellPointers(mesh_p);

  int element_idx = 0;
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
    const lf::geometry::Geometry *geo_p = cell->Geometry();
    const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

    Eigen::Vector2d cur_midpoint = AdvectionFV2D::barycenter(corners);

    std::cout << "Element NR: " << element_idx << std::endl;
    std::cout << "Corners: \n" << corners << std::endl;
    std::cout << "Barycenter: \n" << cur_midpoint << std::endl;
    std::cout << "Normal Vectors\n" << (*normal_vectors)(*cell) << std::endl;
    std::cout << "Area\n" << lf::geometry::Volume(*geo_p) << std::endl;

    for (auto cur_edge : cell->SubEntities(1)){
      const lf::geometry::Geometry* geo_p_edge = cur_edge->Geometry();
      double edge_length = lf::geometry::Volume(*geo_p_edge);
      std::cout << "Edge\n" << lf::geometry::Volume(*geo_p_edge) << std::endl;
    }

    for (const lf::mesh::Entity *neighbour_cell: (*adjacentCells)(*cell)){
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
  const lf::assemble::UniformFEDofHandler dofh(
      mesh_p, {{lf::base::RefEl::kPoint(), 0},
               {lf::base::RefEl::kSegment(), 0},
               {lf::base::RefEl::kTria(), 1},
               {lf::base::RefEl::kQuad(), 1}});

  Eigen::SparseMatrix<double> B_Matrix =
      AdvectionFV2D::initializeMOLODEMatrix(dofh, beta);

  std::cout << "B_Matrix\n" << B_Matrix << std::endl;
  */

  // Testing Output minH
  /*
  double minH = AdvectionFV2D::computeHmin(mesh_p);
  std::cout << "minh " << minH << std::endl;
  */

  // Test of Task 8-8.n
  // Eigen::VectorXd result = AdvectionFV2D::simulateAdvection(dofh);

  // Task 8-8.o
  // Generate a mesh hierarchy
  double T_max = 1.0;

  auto mesh_seq_p{
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, 5)};

  int num_meshes = mesh_seq_p->NumLevels();

  std::vector<int> array_num_cells;
  std::vector<double> array_l1error;
  std::vector<double> array_l2error;
  std::vector<double> array_hmin_inv;
  std::vector<double> array_volume_ref;
  std::vector<double> array_volume_sol;

  // Iterate over all mesh levels
  for (int level = 0; level < num_meshes; level++) {
    std::cout << "Computing L2Error for level: " << level << std::endl;

    // Get the current mesh
    auto cur_mesh = mesh_seq_p->getMesh(level);

    // Create a DOF Hander for the current mesh
    const lf::assemble::UniformFEDofHandler cur_dofh(
        cur_mesh, {{lf::base::RefEl::kPoint(), 0},
                   {lf::base::RefEl::kSegment(), 0},
                   {lf::base::RefEl::kTria(), 1},
                   {lf::base::RefEl::kQuad(), 1}});

    // Get Result from Simulation
    Eigen::VectorXd result = AdvectionFV2D::simulateAdvection(cur_dofh);

    // Get Exact Result at barycenters of cells
    Eigen::VectorXd ref_solution = AdvectionFV2D::refSolution(cur_dofh);

    // Writing vtk files
    std::string sol = "sol";
    std::string ref = "ref";
    std::string f_end = ".vtk";
    std::ostringstream sol_st;
    std::ostringstream ref_st;
    std::ostringstream sol_st_file;
    std::ostringstream ref_st_file;
    sol_st << sol << level;
    ref_st << ref << level;
    sol_st_file << sol << level << f_end;
    ref_st_file << ref << level << f_end;

    lf::io::VtkWriter vtk_writer1(cur_dofh.Mesh(), sol_st_file.str());
    auto cell_data_ref1 =
        lf::mesh::utils::make_CodimMeshDataSet<double>(cur_dofh.Mesh(), 0);
    for (const lf::mesh::Entity *cell : cur_dofh.Mesh()->Entities(0)) {
      int row = cur_dofh.GlobalDofIndices(*cell)[0];
      cell_data_ref1->operator()(*cell) = result[row];
    }
    vtk_writer1.WriteCellData(sol_st.str(), *cell_data_ref1);

    lf::io::VtkWriter vtk_writer2(cur_dofh.Mesh(), ref_st_file.str());
    auto cell_data_ref2 =
        lf::mesh::utils::make_CodimMeshDataSet<double>(cur_dofh.Mesh(), 0);
    for (const lf::mesh::Entity *cell : cur_dofh.Mesh()->Entities(0)) {
      int row = cur_dofh.GlobalDofIndices(*cell)[0];
      cell_data_ref2->operator()(*cell) = ref_solution[row];
    }
    vtk_writer2.WriteCellData(ref_st.str(), *cell_data_ref2);

    // Compute L2 Error in barycenter (bad results)
    /*
    double l2_error = 0;
    for (const lf::mesh::Entity *cell : cur_mesh->Entities(0)) {
      const lf::geometry::Geometry *geo_p = cell->Geometry();
      double area = lf::geometry::Volume(*geo_p);
      int idx = cur_dofh.GlobalDofIndices(*cell)[0];
      l2_error += std::pow((result[idx] - ref_solution[idx]), 2) * area;
    }
    */

    // Compute L1 and L2 Error using high order quadrature
    Eigen::Vector2d x0(0.8, 0.2);
    double d = 0.2;
    auto u0 = [x0, d](Eigen::Vector2d x) -> double {
      double dist = (x - x0).norm();
      if (dist < d) {
        return std::pow(std::cos(M_PI / (2.0 * d) * dist), 2);
      } else {
        return 0.0;
      }
    };
    Eigen::Matrix2d phi_inv;
    double t_mod = T_max / std::sqrt(2.0);
    phi_inv << std::cos(t_mod), std::sin(t_mod), -std::sin(t_mod),
        std::cos(t_mod);

    int quadr_order = 9;
    // Quadrature for reference triangle
    const auto ref_triangle = lf::base::RefEl(lf::base::RefElType::kTria);
    const auto quad_rule_kTria =
        lf::quad::make_QuadRule(ref_triangle, quadr_order);

    // Quadrature for reference cube
    const auto ref_cube = lf::base::RefEl(lf::base::RefElType::kQuad);
    const auto quad_rule_kQuad = lf::quad::make_QuadRule(ref_cube, quadr_order);

    double l1_error = 0;
    double l2_error = 0;
    double volume_ref = 0;
    double volume_sol = 0;
    for (const lf::mesh::Entity *cell : cur_mesh->Entities(0)) {
      const lf::geometry::Geometry *geo_p = cell->Geometry();
      double area = lf::geometry::Volume(*geo_p);
      int idx = cur_dofh.GlobalDofIndices(*cell)[0];

      const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

      if (corners.cols() == 3) {
        const int P = quad_rule_kTria.NumPoints();
        const Eigen::MatrixXd zeta_ref{quad_rule_kTria.Points()};
        const Eigen::MatrixXd zeta{geo_p->Global(zeta_ref)};
        const Eigen::VectorXd w_ref{quad_rule_kTria.Weights()};
        const Eigen::VectorXd gram_dets{geo_p->IntegrationElement(zeta_ref)};

        for (int l = 0; l < P; ++l) {
          l1_error += w_ref[l] *
                      std::abs((u0(phi_inv * zeta.col(l)) - result[idx])) *
                      gram_dets[l];
          l2_error += w_ref[l] *
                      std::pow((u0(phi_inv * zeta.col(l)) - result[idx]), 2) *
                      gram_dets[l];
          volume_ref +=
              w_ref[l] * std::abs(u0(phi_inv * zeta.col(l))) * gram_dets[l];
          volume_sol += w_ref[l] * std::abs(result[idx]) * gram_dets[l];
        }
      } else if (corners.cols() == 4) {
        const int P = quad_rule_kQuad.NumPoints();
        const Eigen::MatrixXd zeta_ref{quad_rule_kQuad.Points()};
        const Eigen::MatrixXd zeta{geo_p->Global(zeta_ref)};
        const Eigen::VectorXd w_ref{quad_rule_kQuad.Weights()};
        const Eigen::VectorXd gram_dets{geo_p->IntegrationElement(zeta_ref)};

        for (int l = 0; l < P; ++l) {
          l1_error += w_ref[l] *
                      std::abs((u0(phi_inv * zeta.col(l)) - result[idx])) *
                      gram_dets[l];
          l2_error += w_ref[l] *
                      std::pow((u0(phi_inv * zeta.col(l)) - result[idx]), 2) *
                      gram_dets[l];
          volume_ref +=
              w_ref[l] * std::abs(u0(phi_inv * zeta.col(l))) * gram_dets[l];
          volume_sol += w_ref[l] * std::abs(result[idx]) * gram_dets[l];
        }
      } else {
        throw std::runtime_error("Error in L2 Geometrie");
      }
    }
    l2_error = std::sqrt(l2_error);

    std::cout << "L2Error at level " << level << " :" << l2_error << std::endl;

    // Store the results in the vectors
    array_num_cells.push_back(cur_dofh.NumDofs());
    array_l1error.push_back(l1_error);
    array_l2error.push_back(l2_error);
    array_volume_sol.push_back(volume_sol);
    array_volume_ref.push_back(volume_ref);

    double hmin_inv = 1. / AdvectionFV2D::computeHmin(cur_mesh);
    array_hmin_inv.push_back(hmin_inv);
  }

  // Task 8-8.q
  std::vector<int> array_clf_thres;

  // Iterate over all mesh levels
  for (int level = 0; level < num_meshes; level++) {
    // Get the current mesh
    auto cur_mesh = mesh_seq_p->getMesh(level);

    // Create a DOF Hander for the current mesh
    const lf::assemble::UniformFEDofHandler cur_dofh(
        cur_mesh, {{lf::base::RefEl::kPoint(), 0},
                   {lf::base::RefEl::kSegment(), 0},
                   {lf::base::RefEl::kTria(), 1},
                   {lf::base::RefEl::kQuad(), 1}});

    std::cout << "Computing threshold for level: " << level << std::endl;

    // Save the CFL threshold to array
    // For debugging of convergence:
    // int threshold = 1;
    int threshold = AdvectionFV2D::findCFLthreshold(cur_dofh, beta, T_max);
    array_clf_thres.push_back(threshold);
  }

  // Write Output file of Task 8-8.o and 8-8.q
  std::ofstream csv_file;
  csv_file.open("advectionfv2d.csv");
  for (int level = 0; level < num_meshes; level++) {
    std::cout << "Cells: " << array_num_cells.at(level)
              << " | L2Error: " << array_l2error.at(level)
              << " | L1Error: " << array_l1error.at(level)
              << " | Vol. Ref: " << array_volume_ref.at(level)
              << " | Vol. Sol: " << array_volume_sol.at(level)
              << " | Thres: " << array_clf_thres.at(level)
              << " | 1/Hmin: " << array_hmin_inv.at(level) << std::endl;
    csv_file << array_num_cells.at(level) << "," << array_hmin_inv.at(level)
             << "," << array_l1error.at(level) << "," << array_l2error.at(level)
             << "," << array_clf_thres.at(level) << "\n";
  }
  csv_file.close();

  // const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
  //    Eigen::DontAlignCols, ", ", "\n");
  // std::cout << v.transpose().format(CSVFormat) << std::endl;
  return 0;
}
