/**
 * @file preptrimesh.h
 * @brief Save p.w. linear function on triangular mesh for visualization with matlab
 * @author R. Hiptmair
 * @date April 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/uscalfe/uscalfe.h>

namespace PrepTriMesh {
  /** @brief Save p.w. linear function on triangular mesh for visualization with matlab
   *
   * @param fes_p pointer to lowest-order Lagrangian FE space
   * @param mu vector of nodal values of p.w. linear function
   * @param filename name of output file (.m appended, if not there already)
   *
   * This generates a file containing a MATLAB function
   *
   * function [x,y,TRI,QUAD,EDS] = filename();
   *
   * calling lf::io::writeMatlab. It creates a second file "nodvals_<filename>.m",
   * containing a vector 'nodevals' of nodal values. 
   */
void prepTriMesh(std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
		 const Eigen::VectorXd& mu, std::string filename);
}
