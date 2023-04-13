/**
 * @file minimalgraphsurface.cc
 * @brief NPDE homework 5-3 Minimal Graph Surface code
 * @author R. Hiptmair
 * @date April 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "preptrimesh.h"

#include <lf/assemble/assembly_types.h>
#include <lf/base/base.h>
#include <lf/base/lf_assert.h>
#include <lf/base/ref_el.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/io/write_matlab.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>

#include <fstream>

namespace PrepTriMesh {
// The only function
void prepTriMesh(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
    const Eigen::VectorXd& mu, std::string filename) {
  // Obtain reference to the underlying finite element mesh
  const lf::mesh::Mesh& mesh{*fes_p->Mesh()};
  // The local-to-global index mapping
  const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
  // Get the number of degrees of freedom = dimension of FE space
  const lf::base::size_type N_dofs(dofh.NumDofs());
  LF_VERIFY_MSG(mesh.DimMesh() == 2, "Only 2D meshes supported");
  LF_VERIFY_MSG(mesh.DimWorld() == 2, "Only planar meshes supported");
  LF_ASSERT_MSG(N_dofs == mesh.NumEntities(2), "No d.o.f.s !+ no. nodes!");
  LF_VERIFY_MSG(mu.size() == N_dofs, "Vector length mismatch!");
  LF_VERIFY_MSG(mesh.NumEntities(lf::base::RefEl::kQuad()) == 0,
                "Only for triangular meshes!");
  // Generate a MATLAB function that provides information about the mesh
  lf::io::writeMatlab(mesh, filename);
  // Create a second .m file containing the nodal values
  std::vector<double> nodvals(N_dofs);
  for (const lf::mesh::Entity* node_ptr : mesh.Entities(2)) {
    nonstd::span<const lf::assemble::gdof_idx_t> dof_idx{
        dofh.GlobalDofIndices(*node_ptr)};
    LF_VERIFY_MSG(dof_idx.size() == 1, "Only 1 dof per node allowed!");
    // Obtain index of current node
    const lf::base::size_type node_idx = mesh.Index(*node_ptr);
    LF_VERIFY_MSG(node_idx < N_dofs, "Node index >= N_dofs!");
    // Obtain value for current node
    double nodeval = mu[dof_idx[0]];
    // Store value for current node
    nodvals[node_idx] = nodeval;
  }
  // Write nodal values to file
  // Preprocess filename: append .m, if not there already
  {
    lf::base::size_type fn_len = filename.size();
    if ((fn_len > 1) && (filename[fn_len - 2] != '.') &&
        (filename[fn_len - 1] != 'm')) {
      filename += ".m";
    }
  }
  // Put 'nodvals' in front of the filename
  filename = "nodvals_" + filename;

  std::ofstream file(filename);
  if (file.is_open()) {
    file << "% Vector of nodal values" << std::endl;
    file << "nodevals = [ ";
    for (double v : nodvals) {
      file << v << " ";
    }
    file << "];" << std::endl;
  }
}
  
}  // namespace PrepTriMesh
