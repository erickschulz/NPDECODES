/**
 * @file
 * @brief NPDE homework SymplecticTimesteppingWaves
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "symplectictimesteppingwaves.h"

namespace SymplecticTimesteppingWaves {

/* SAM_LISTING_BEGIN_7 */
void wavePropSimulation(unsigned int m) {
  //====================
  // Your code goes here
  //====================
}

/* SAM_LISTING_END_7 */

void progress_bar::write(double fraction) {
  // clamp fraction to valid range [0,1]
  if (fraction < 0)
    fraction = 0;
  else if (fraction > 1)
    fraction = 1;

  auto width = bar_width - message.size();
  auto offset = bar_width - static_cast<unsigned>(width * fraction);

  os << '\r' << message;
  os.write(full_bar.data() + offset, width);
  os << " [" << std::setw(3) << static_cast<int>(100 * fraction) << "%] "
     << std::flush;
}

/* SAM_LISTING_BEGIN_5 */
double testStab() {
  double maxUniformTimestep;
  //====================
  // Your code goes here
  //====================
  return maxUniformTimestep;
}
/* SAM_LISTING_END_5 */

}  // namespace SymplecticTimesteppingWaves
