/**
 * @file transpsemilagr_main.cc
 * @brief NPDE homework TranspSemiLagr Main file
 * @author Philippe Peter
 * @date November 2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <filesystem>
#include <memory>

#include "transpsemilagr.h"

int main() {
  // The equation is solved on the test mesh circle.msh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory),
                            CURRENT_SOURCE_DIR "/../meshes/circle.msh");
  auto mesh_p = reader.mesh();

  // construct linear finite element space
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // initial conditions:
  // const auto u0 = [](const Eigen::Vector2d& x){ return 1-(x(0)* x(0) +
  // x(1)*x(1));};
  const auto u0 = [](const Eigen::Vector2d& x) {
    if (x(1) * x(1) + x(0) * x(0) > 0.99) {
      return 0.0;
    } else {
      return -x(1) - x(0);
    }
  };

  Eigen::VectorXd u0_vector = lf::fe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u0));

  // compute the solutions after 1st timestep
  Eigen::VectorXd sol_rot_1 =
      TranspSemiLagr::solverot(fe_space, u0_vector, 1, 0.1);
  lf::fe::MeshFunctionFE mf_sol_rot_1(fe_space, sol_rot_1);
  Eigen::VectorXd sol_trp_1 =
      TranspSemiLagr::solvetrp(fe_space, u0_vector, 1, 0.1);
  lf::fe::MeshFunctionFE mf_sol_trp_1(fe_space, sol_trp_1);

  // compute solutions at final time
  Eigen::VectorXd sol_rot_10 =
      TranspSemiLagr::solverot(fe_space, u0_vector, 10, 1.0);
  lf::fe::MeshFunctionFE mf_sol_rot_10(fe_space, sol_rot_10);
  Eigen::VectorXd sol_trp_10 =
      TranspSemiLagr::solvetrp(fe_space, u0_vector, 10, 1.0);
  lf::fe::MeshFunctionFE mf_sol_trp_10(fe_space, sol_trp_10);

  // OUTPUT RESULTS
  // create directory `transp_semi_lagr_solution` if it doesn't exist yet:
  if (!std::filesystem::is_directory("transp_semi_lagr_solution") ||
      !std::filesystem::exists("transp_semi_lagr_solution")) {
    std::filesystem::create_directory("transp_semi_lagr_solution");
  }
  // construct writers
  lf::io::VtkWriter vtk_writer_rot_1(
      mesh_p, CURRENT_BINARY_DIR "/transp_semi_lagr_solution/rot_1.vtk");
  lf::io::VtkWriter vtk_writer_rot_10(
      mesh_p, CURRENT_BINARY_DIR "/transp_semi_lagr_solution/rot_10.vtk");
  lf::io::VtkWriter vtk_writer_trp_1(
      mesh_p, CURRENT_BINARY_DIR "/transp_semi_lagr_solution/trp_1.vtk");
  lf::io::VtkWriter vtk_writer_trp_10(
      mesh_p, CURRENT_BINARY_DIR "/transp_semi_lagr_solution/trp_10.vtk");

  // output data
  vtk_writer_rot_1.WritePointData("rot_1", mf_sol_rot_1);
  vtk_writer_rot_10.WritePointData("rot_10", mf_sol_rot_10);
  vtk_writer_trp_1.WritePointData("trp_1", mf_sol_trp_1);
  vtk_writer_trp_10.WritePointData("trp_10", mf_sol_trp_10);

  return 0;
}
