/**
 * @ file
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Core>
#include <iostream>

#include "debuggingwithlldb.h"

/* SAM_LISTING_BEGIN_1 */
int main(int /*argc*/, char** /*argv*/) {
  // A function using some of LehrFEM++'s facilities
  DebuggingWithLLDB::ReadAndOutputMesh("../ljoint.msh");
  return 0;
}
/* SAM_LISTING_END_1 */
