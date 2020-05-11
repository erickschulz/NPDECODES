#ifndef __LINFEREACTDIFF_H
#define __LINFEREACTDIFF_H
/**
 This homework problem consists of reading a simple, gmesh generated, mesh on
 the unit square and solving a simple reaction diffusion system using LehrFEM++
 */

#include <lf/base/base.h>
#include <lf/mesh/mesh.h>
#include <lf/refinement/refinement.h>

#include <Eigen/Core>
#include <memory>

namespace LinFeReactDiff {

std::shared_ptr<lf::refinement::MeshHierarchy> generateMeshHierarchy(
    const lf::base::size_type levels);

Eigen::VectorXd solveFE(std::shared_ptr<const lf::mesh::Mesh> mesh);

double computeEnergy(std::shared_ptr<const lf::mesh::Mesh> mesh,
                     Eigen::VectorXd mu);
}  // namespace LinFeReactDiff
#endif  // define __LINFEREACTDIFF_H
