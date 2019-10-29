/*
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss
 * @date   13.03.2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <boost/filesystem.hpp>
#include <iomanip>
#include "l2_error_cr_discretization_dirichlet_bvp.h"

using namespace NonConformingCrouzeixRaviartFiniteElements;

int main() {
  // Obtain paths to mesh files
  boost::filesystem::path file_path = __FILE__;
  auto mesh_dir_path = file_path.parent_path().parent_path() / "meshes";

  auto sep = std::setw(15);
  std::cout << std::left << sep << "Number of dofs"
            << " | "
            << "L2 error" << std::endl;

  // Loop over meshes and output number of dofs and L2 norm of error
  for (int i = 1; i <= 4; ++i) {
    std::string mesh_path =
        (mesh_dir_path / ("refined_square" + std::to_string(i) + ".msh"))
            .string();

    // Read mesh from file
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    const lf::io::GmshReader reader(std::move(mesh_factory), mesh_path);
    auto mesh_ptr = reader.mesh();

    // TODO: 2-14.i)
    /* BEGIN_SOLUTION */
    /* TODO Your implementation goes here! */
    /* END_SOLUTION */

    // TODO: 2-14.y)
    /* BEGIN_SOLUTION */
    /* TODO Your implementation goes here! */
    /* END_SOLUTION */
  }

  return 0;
}
