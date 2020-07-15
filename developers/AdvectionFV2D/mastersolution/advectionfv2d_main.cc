/**
 * @file advectionfv2d_main.cc
 * @brief NPDE homework AdvectionFV2D code
 * @author Philipp Egg
 * @date 21.06.2020
 * @copyright Developed at ETH Zurich
 */

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>

#include "advectionfv2d.h"

/* SAM_LISTING_BEGIN_1 */
int main() {
#if SOLUTION
  // Define velocity field beta
  // Note that the problem description requires ||B|| <= 1
  auto beta = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(-x[1], x[0]) / std::sqrt(2.0);
  };

  // Functor for initial bump
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

  // Task 8-8.o
  // Generate a mesh hierarchy
  double T = 1.0;

  //////////////////////////////////////////////////////////////////////////////
  // TODO inconsistancy: g vs. G
  //////////////////////////////////////////////////////////////////////////////
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1. / 3.);

  auto mesh_seq_p{
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, 6)};

  std::vector<int> vector_num_cells;
  std::vector<double> vector_l2error;

  // Iterate over mesh levels starting from thrid refinement
  int num_meshes = mesh_seq_p->NumLevels();
  for (int level = 3; level < num_meshes; ++level) {
    std::cout << "Computing L2Error for level: " << level << std::endl;

    // Get the current mesh
    auto cur_mesh = mesh_seq_p->getMesh(level);

    // Create a DOF Hander for the current mesh
    const lf::assemble::UniformFEDofHandler cur_dofh(
        cur_mesh, {{lf::base::RefEl::kPoint(), 0},
                   {lf::base::RefEl::kSegment(), 0},
                   {lf::base::RefEl::kTria(), 1},
                   {lf::base::RefEl::kQuad(), 1}});

    // Compute cell normals
    std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
        Eigen::Matrix<double, 2, Eigen::Dynamic>>>
        normal_vectors = AdvectionFV2D::computeCellNormals(cur_dofh.Mesh());

    // Compute adjecent cells
    std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
        std::array<const lf::mesh::Entity *, 4>>>
        adjacentCells = AdvectionFV2D::getAdjacentCellPointers(cur_dofh.Mesh());

    // Get result from simulation
    Eigen::VectorXd result = AdvectionFV2D::simulateAdvection(
        cur_dofh, beta, u0, adjacentCells, normal_vectors, T);

    // Get exact result at barycenters of cells
    Eigen::VectorXd ref_solution = AdvectionFV2D::refSolution(cur_dofh, u0, T);

    // Compute L2 error in barycenter
    double l2_error = 0;
    for (const lf::mesh::Entity *cell : cur_mesh->Entities(0)) {
      const lf::geometry::Geometry *geo_p = cell->Geometry();
      double area = lf::geometry::Volume(*geo_p);
      int idx = cur_dofh.GlobalDofIndices(*cell)[0];
      l2_error += std::pow((result[idx] - ref_solution[idx]), 2) * area;
    }
    l2_error = std::sqrt(l2_error);
    vector_num_cells.push_back(cur_dofh.NumDofs());
    vector_l2error.push_back(l2_error);
    std::cout << "L2Error at level " << level << ": " << l2_error << std::endl;

    // Writing vtk files (optional part)
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
    auto cell_data_ref =
        lf::mesh::utils::make_CodimMeshDataSet<double>(cur_dofh.Mesh(), 0);
    for (const lf::mesh::Entity *cell : cur_dofh.Mesh()->Entities(0)) {
      int row = cur_dofh.GlobalDofIndices(*cell)[0];
      cell_data_ref->operator()(*cell) = result[row];
    }
    vtk_writer1.WriteCellData(sol_st.str(), *cell_data_ref);

    lf::io::VtkWriter vtk_writer2(cur_dofh.Mesh(), ref_st_file.str());
    auto cell_data_sol =
        lf::mesh::utils::make_CodimMeshDataSet<double>(cur_dofh.Mesh(), 0);
    for (const lf::mesh::Entity *cell : cur_dofh.Mesh()->Entities(0)) {
      int row = cur_dofh.GlobalDofIndices(*cell)[0];
      cell_data_sol->operator()(*cell) = ref_solution[row];
    }
    vtk_writer2.WriteCellData(ref_st.str(), *cell_data_sol);
  }

  // Task 8-8.q
  // Compute threshold for fourth refinement level
  int level = 4;
  auto cur_mesh = mesh_seq_p->getMesh(level);

  // Create a DOF Hander for the current mesh
  const lf::assemble::UniformFEDofHandler cur_dofh(
      cur_mesh, {{lf::base::RefEl::kPoint(), 0},
                 {lf::base::RefEl::kSegment(), 0},
                 {lf::base::RefEl::kTria(), 1},
                 {lf::base::RefEl::kQuad(), 1}});

  int threshold = AdvectionFV2D::findCFLthreshold(cur_dofh, beta, T);
  int cfl_thres = int((T / AdvectionFV2D::computeHmin(cur_mesh) + 1));
  std::cout << "Threshold for level: " << level << " is: " << threshold
            << " | Threshold from CFL is: " << cfl_thres << std::endl;

  // Write Output file of Task 8-8.o and 8-8.q
  std::ofstream csv_file;
  csv_file.open("advectionfv2d.csv");
  for (int i = 0; i < vector_num_cells.size(); ++i) {
    std::cout << "Cells: " << vector_num_cells.at(i)
              << " | L2Error: " << vector_l2error.at(i) << std::endl;
    csv_file << vector_num_cells.at(i) << "," << vector_l2error.at(i) << "\n";
  }
  csv_file.close();

  // Print convergence rates
  for (int i = 1; i < vector_num_cells.size(); ++i) {
    double conv_rate =
        (std::log(vector_l2error[i - 1]) - std::log(vector_l2error[i])) /
        (std::log(vector_num_cells[i]) - std::log(vector_num_cells[i - 1]));
    std::cout << "Conv. Rate " << i << "-" << i + 1 << "  is: " << conv_rate
              << std::endl;
  }

  // Plot results from 8-8.o
  std::system("python3 " CURRENT_SOURCE_DIR
              "/advectionfv2d.py " CURRENT_BINARY_DIR
              "/advectionfv2d.csv " CURRENT_BINARY_DIR "/solution.eps");

  return 0;
#else
  //====================
  // Your code goes here
  //====================
  return 0;
#endif
}
/* SAM_LISTING_END_1 */
