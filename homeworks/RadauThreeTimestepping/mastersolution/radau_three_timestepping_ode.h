/** @file
 * @brief NPDE RadauThreeTimestepping
 * @author Erick Schulz
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include <math.h>
#include <iostream>
#include <vector> 

namespace RadauThreeTimestepping {

std::vector<double> twoStageRadauTimesSteppingLinScalODE(unsigned int);
void testConvergenceTwoStageRadauLinScalODE();

}  // RadauThreeTimestepping
