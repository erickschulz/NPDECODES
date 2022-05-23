/**
 * @file upwindfinitevolume_main.cc
 * @brief NPDE homework UpwindFiniteVolume code
 * @author Philipp Egg
 * @date 08.09.2020
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/gmsh_reader.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/refinement/refinement.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "upwindfinitevolume.h"

int main() {
  /* SAM_LISTING_BEGIN_1 */
  // Read in mesh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory),
                            CURRENT_SOURCE_DIR "/mesh.msh");

  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader.mesh();

  auto mesh_seq_p{
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, 6)};
  int num_meshes = mesh_seq_p->NumLevels();

  //====================
  // Your code goes here
  //====================
  /* SAM_LISTING_END_1 */
  return 0;
}
