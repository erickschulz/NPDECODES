/**
 * @ file boundarylength.cc
 * @ brief NPDE homework LengthOfBoundary code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "boundarylength.h"

#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>

namespace LengthOfBoundary {

/* SAM_LISTING_BEGIN_1 */
double volumeOfDomain(const std::shared_ptr<lf::mesh::Mesh> mesh_p) {
  double volume = 0.0;
  //====================
  // Your code goes here
  //====================

  return volume;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double lengthOfBoundary(const std::shared_ptr<lf::mesh::Mesh> mesh_p) {
  double length = 0.0;
  //====================
  // Your code goes here
  //====================

  return length;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::pair<double, double> measureDomain(std::string filename) {
  double volume, length;

  //====================
  // Your code goes here
  //====================

  return {volume, length};
}
/* SAM_LISTING_END_3 */

}  // namespace LengthOfBoundary
