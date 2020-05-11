/*
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss, edited Amélie Loher
 * @date   13.03.2019, 03.03.20
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>

#include <iomanip>

#include "crl2errordirichletbvp.h"

using namespace NonConformingCrouzeixRaviartFiniteElements;

int main() {
  auto sep = std::setw(15);
  std::cout << std::left << sep << "Number of dofs"
            << " | "
            << "L2 error" << std::endl;

  // Loop over meshes and output number of dofs and L2 norm of error
  for (int i = 1; i <= 4; ++i) {
    std::string mesh_file = CURRENT_SOURCE_DIR "/../meshes/refined_square" +
                            std::to_string(i) + ".msh";

    // Read mesh from file
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
    auto mesh_p = reader.mesh();

// TODO: task 2-14.h)
#if SOLUTION
    // Initialize CR-DofHandler from mesh
    lf::assemble::UniformFEDofHandler dof_handler(
        mesh_p, {{lf::base::RefEl::kPoint(), 0},
                 {lf::base::RefEl::kSegment(), 1},
                 {lf::base::RefEl::kTria(), 0},
                 {lf::base::RefEl::kQuad(), 0}});
#else
    //====================
    // Your code goes here
    //====================
#endif

// TODO: task 2-14.y)
#if SOLUTION
    std::cout << std::left << sep << dof_handler.NumDofs() << " | "
              << L2errorCRDiscretizationDirichletBVP(mesh_file) << std::endl;
#else
    //====================
    // Your code goes here
    //====================
#endif
  }

  return 0;
}
