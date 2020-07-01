/**
 * @file qfeinterpolator.h
 * @brief NPDE homework DebuggingFEM code
 * @author Simon Meierhans
 * @date 27/03/2019
 * @copyright Developed at ETH Zurich
 */

#ifndef NPDECODES_DEBUGGINGFEM_QFEINTERPOLATOR_H_
#define NPDECODES_DEBUGGINGFEM_QFEINTERPOLATOR_H_

#include <Eigen/Core>

#include <lf/assemble/assemble.h>
#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

namespace DebuggingFEM {

using size_type = lf::base::size_type;

/**
 * @brief get the global coordinates of a interpolation point in a cell
 * @param idx local index
 *        cell the cell to be used
 * @returns The global coordinate of the i-th interpolation node in the given
 * cell
 */
Eigen::Vector2d globalCoordinate(int idx, const lf::mesh::Entity &cell);

/**
 * @brief interpolate function over a second order lagrangian finite element
 * space
 * @param dofh dof-handler
 *        f function to interpolate
 * @returns A vector containing thebasis function coefficients
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR>
Eigen::VectorXd interpolateOntoQuadFE(const lf::assemble::DofHandler &dofh,
                                      FUNCTOR &&f) {
  // get mesh and set up finite element space
  auto mesh = dofh.Mesh();

  // initialize result vector
  const size_type N_dofs(dofh.NumDofs());
  Eigen::VectorXd result = Eigen::VectorXd::Zero(N_dofs);

#if SOLUTION
  for (const lf::mesh::Entity *cell : mesh->Entities(0)) {
    // get local to global map for the cell
    auto glob_ind = dofh.GlobalDofIndices(*cell);
    // This is a rather clumsy implementation. A much better way to accomplish
    // this is to pass a full matrix of reference coorindates to the Global()
    // methof of the current cell entity. This will immediately give the global
    // coorindates of all local interpolation nodes and we can dispense with the
    // auxiliary function. Afterwards we can loop over the columns of the
    // matrix.
    for (int i = 0; i < 6; i++) {
      // update the result vector
      auto coords = globalCoordinate(i, *cell);
      result(glob_ind[i]) = f(coords);
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return result;
}
/* SAM_LISTING_END_1 */

} // namespace DebuggingFEM

#endif
