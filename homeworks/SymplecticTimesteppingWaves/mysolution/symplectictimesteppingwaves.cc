/**
 * @file
 * @brief NPDE homework SymplecticTimesteppingWaves
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "symplectictimesteppingwaves.h"

namespace SymplecticTimesteppingWaves {

void wavePropSimulation(unsigned int m) {
  /* SOLUTION_BEGIN */
  /* TODO Your implementation goes here! */
  /* SOLUTION_END */
}

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

double testStab(){
  double maxUniformTimestep;
  /* SOLUTION_BEGIN */
  /* TODO Your implementation goes here! */
  /* SOLUTION_END */
  return maxUniformTimestep;
}

}  // namespace SymplecticTimesteppingWaves
