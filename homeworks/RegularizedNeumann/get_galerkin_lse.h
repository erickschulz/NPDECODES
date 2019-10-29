#ifndef __REGNEUMANN_H
#define __REGNEUMANN_H
/**
 * @ file regNeumann.h
 * @ brief NPDE homework RegularizedNeumann code
 * @ author Christian Mitsch
 * @ date 11.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include <iostream>
#include <memory>

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace RegularizedNeumann {

/* SAM_LISTING_BEGIN_1 */
template <typename FUNCT_F, typename FUNCT_H>
std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> getGalerkinLSE(
    const std::shared_ptr<lf::uscalfe::ScalarUniformFESpace<double>> fe_space,
    const FUNCT_F &f, const FUNCT_H &h) {
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  const std::size_t N_dofs = dofh.NoDofs();
  // Right-hand-side vector; don't forget to set to zero initially!
  Eigen::VectorXd rhs_vec(N_dofs);
  rhs_vec.setZero();

  // ASSEMBLE MATRIX
  lf::assemble::COOMatrix<double> A_aux(N_dofs, N_dofs);
  lf::uscalfe::LinearFELaplaceElementMatrix my_mat_provider{};
  // co-dimension 0 because we locally assemble on cells
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, my_mat_provider, A_aux);

  // ASSEMBLE GLOBAL RHS VECOR
  // Volume part $v\mapsto\int_{\Omega}fv\,\mathrm{d}\Bx$ of the rhs linear
  // functional
  lf::uscalfe::ScalarLoadElementVectorProvider my_vec_provider(fe_space, f);
  // co-dimension 0 because we locally assemble on cells
  lf::assemble::AssembleVectorLocally(0, dofh, my_vec_provider, rhs_vec);

  // Boundary part $v\mapsto\int_{\partial\Omega}hv\,\mathrm{d}S$ of the rhs
  // linear functional. Only edges on the boundary should be included in
  // assembly, which is achieved by passing a suitable selector predicate to the
  // builder object through a lambda function
  auto bd_edges{lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 1)};
  lf::uscalfe::ScalarLoadEdgeVectorProvider my_vec_provider_edge(
      fe_space, h, [&bd_edges](const lf::mesh::Entity &edge) -> bool {
        return bd_edges(edge);
      });
  // co-dimension 1 because we locally assemble on edges !
  lf::assemble::AssembleVectorLocally(1, dofh, my_vec_provider_edge, rhs_vec);

  return std::make_pair(A_aux.makeSparse(), rhs_vec);
}
/* SAM_LISTING_END_1 */

}  // namespace RegularizedNeumann

#endif
