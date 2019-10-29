/**
 * @ file boundarylength.cc
 * @ brief NPDE homework LengthOfBoundary code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "boundarylength.h"

namespace LengthOfBoundary {

/* SAM_LISTING_BEGIN_1 */
double volumeOfDomain(const std::shared_ptr<lf::mesh::Mesh> mesh) {
  double volume = 0.0;
  /* BEGIN_SOLUTION */
  // iterate over all cells (co-dimension = 0)
  for (const lf::mesh::Entity &ent : mesh->Entities(0)) {
    lf::geometry::Geometry *geo_ptr = ent.Geometry();
    volume += lf::geometry::Volume(*geo_ptr);
  }
  /* END_SOLUTION */
  return volume;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double lengthOfBoundary(const std::shared_ptr<lf::mesh::Mesh> mesh) {
  double length = 0.0;
  /* BEGIN_SOLUTION */
  // This function returns an array of flags
  // a flag is true if an edge is part of the boundary
  auto edge_marker = lf::mesh::utils::flagEntitiesOnBoundary(mesh, 1);

  // iterate over all edges (co-dimension = 1)
  for (const lf::mesh::Entity &ent : mesh->Entities(1)) {
    // check if edge is part of the boundary
    if (edge_marker(ent)) {
      lf::geometry::Geometry *geo_ptr = ent.Geometry();
      length += lf::geometry::Volume(*geo_ptr);
    }
  }
  /* END_SOLUTION */
  return length;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::pair<double, double> measureDomain(std::string msh_file_name) {
  double volume, length;
  /* BEGIN_SOLUTION */
  // read in mesh file
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), msh_file_name);
  auto mesh_p = reader.mesh();

  // call the functions we already implemented
  length = LengthOfBoundary::lengthOfBoundary(mesh_p);
  volume = LengthOfBoundary::volumeOfDomain(mesh_p);
  /* END_SOLUTION */

  return {volume, length};
}
/* SAM_LISTING_END_3 */

}  // namespace LengthOfBoundary
