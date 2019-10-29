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
  // TODO Your implementation goes here!
  /* END_SOLUTION */
  return volume;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double lengthOfBoundary(const std::shared_ptr<lf::mesh::Mesh> mesh) {
  double length = 0.0;
  /* BEGIN_SOLUTION */
  // TODO Your implementation goes here!
  /* END_SOLUTION */
  return length;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::pair<double, double> measureDomain(std::string msh_file_name) {
  double volume, length;
  /* BEGIN_SOLUTION */
  // TODO Your implementation goes here!
  /* END_SOLUTION */

  return {volume, length};
}
/* SAM_LISTING_END_3 */

}  // namespace LengthOfBoundary
