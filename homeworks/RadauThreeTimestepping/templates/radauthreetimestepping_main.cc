/**
 * @file radauthreetimestepping_main.cc
 * @brief NPDE homework RadauThreeTimestepping
 * @author Erick Schulz
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include "radauthreetimestepping.h"
#include "radauthreetimesteppingode.h"

#include <iostream>
#include <memory>

#include <Eigen/Core>

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>

using namespace RadauThreeTimestepping;

int main(int /*argc*/, char ** /*argv*/) {
  /* Solving the ODE problem */
  // This function prints to the terminal the convergence rates and average rate
  // of a convergence study performed for the ODE (d/dt)y = -y.
  testConvergenceTwoStageRadauLinScalODE();

  /* Solving the parabolic heat equation */
  // Create a Lehrfem++ square tensor product mesh
  lf::mesh::hybrid2d::TPTriagMeshBuilder builder(
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
  // Set mesh parameters following the Builder pattern
  // Domain is the unit square
  builder.setBottomLeftCorner(Eigen::Vector2d{-1.0, -1.0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNumXCells(50)
      .setNumYCells(50);
  auto mesh_p = builder.Build();

  /* SAM_LISTING_BEGIN_1 */
  //====================
  // Your code goes here
  //====================

  return 0;
}
