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
typedef Eigen::Triplet<double> Triplet;

// SUB-EXERCISE c)

template <typename FUNCT_F, typename FUNCT_H>
std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> getGalerkinLSE_dropDof(
    const std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>>
        fe_space,
    const FUNCT_F f, FUNCT_H h) {
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  const std::size_t N_dofs = dofh.NoDofs();
  // Right-hand-side vector; don't forget to set to zero initially!
  Eigen::VectorXd rhs_vec(N_dofs);
  rhs_vec.setZero();
  // For returning the sparse Galerkin matrix
  Eigen::SparseMatrix<double> A(N_dofs, N_dofs);

  // ASSEMBLE MATRIX
  lf::assemble::COOMatrix<double> A_aux(N_dofs, N_dofs);
  lf::uscalfe::LinearFELaplaceElementMatrix my_mat_provider{};
  // co-dimension 0 because we locally assemble on cells
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, my_mat_provider, A_aux);

  // ASSEMBLE GLOBAL RHS VECOR
  // Volume part $v\mapsto\int_{\Omega}fv\,\mathrm{d}\Bx$ of the rhs linear
  // functional
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(f)>
      my_vec_provider(fe_space, f);
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

  // Now fix the solution to be 0 at the node p
  /* BEGIN_SOLUTION */
  // Tell the function fix_flagged_solution_components that is should
  // modify the matrix and the right-hand-side vector in a way that amounts to
  // dropping the first tent function from the basis.
  /* SAM_LISTING_BEGIN_1 */
  auto selector = [](lf::base::glb_idx_t idx) -> std::pair<bool, double> {
    if (idx == 0) {
      return {true, 0.0};  // fix first d.o.f. to zero
    } else {
      return {false, 42.0};  // keep all others
    }
  };
  lf::assemble::fix_flagged_solution_components(selector, A_aux, rhs_vec);
  /* SAM_LISTING_END_1 */
  /* END_SOLUTION */
  A = A_aux.makeSparse();
  return std::make_pair(A, rhs_vec);
}

// SUB-EXERCISE f)

// You can use write a helper class which should implement
// ENTITY_VECTOR_PROVIDER to calculate vector c using LehrFEM assembly functions
/* SAM_LISTING_BEGIN_6 */
class VecHelper {
 public:
  explicit VecHelper() {}
  bool isActive(const lf::mesh::Entity &entity) const { return true; }
  Eigen::Vector3d Eval(const lf::mesh::Entity &entity) {
    LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                  "Function only defined for triangular cells");
    Eigen::Vector3d result;
    /* BEGIN_SOLUTION */
    // Obtain shape information for the cell
    const lf::geometry::Geometry *geo_ptr = entity.Geometry();
    // Fetch area
    const double area = lf::geometry::Volume(*geo_ptr);
    // Initialize element vector |K|/3*[1,1,1]^T
    result = (area / 3.0) * Eigen::Vector3d::Ones();
    /* END_SOLUTION */
    return result;
  }
};
/* SAM_LISTING_END_6 */

/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd assembleVector_c(const lf::assemble::DofHandler &dofh) {
  Eigen::VectorXd c(dofh.NoDofs());
  /* BEGIN_SOLUTION */
  // Do not forget to initialize vector before assembly!
  c.setZero();
  // ELEMENT_VECTOR_BUILDER object
  RegularizedNeumann::VecHelper my_vec_provider_c{};
  // Cell (= codim-0 entities)-oriented assembly into c
  lf::assemble::AssembleVectorLocally(0, dofh, my_vec_provider_c, c);
  /* END_SOLUTION */
  return c;
}
/* SAM_LISTING_END_5 */

template <typename FUNCT_F, typename FUNCT_H>
std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> getGalerkinLSE_augment(
    const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    const FUNCT_F f, const FUNCT_H h) {
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  const std::size_t N_dofs = dofh.NoDofs() + 1;
  Eigen::VectorXd rhs_vec(N_dofs);
  rhs_vec.setZero();
  Eigen::SparseMatrix<double> A(N_dofs, N_dofs);

  // ASSEMBLE MATRIX
  lf::assemble::COOMatrix<double> A_aux(N_dofs, N_dofs);
  lf::uscalfe::LinearFELaplaceElementMatrix my_mat_provider{};
  // co-dimension 0 because we locally assemble on cells
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, my_mat_provider, A_aux);

  // ASSEMBLE GLOBAL RHS VECOR
  // Volume part $v\mapsto\int_{\Omega}fv\,\mathrm{d}\Bx$ of the rhs linear
  // functional
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(f)>
      my_vec_provider(fe_space, f);
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
  // co-dimension 1 because we locally assemble on edges
  lf::assemble::AssembleVectorLocally(1, dofh, my_vec_provider_edge, rhs_vec);

  // Calculate c and insert it into A
  /* BEGIN_SOLUTION */
  /* SAM_LISTING_BEGIN_2 */
  // Calculate vector c
  Eigen::VectorXd c = assembleVector_c(dofh);
  // Add c to the matrix
  for (int i = 0; i < c.size(); i++) {
    // add triplets for vector c
    A_aux.AddToEntry(N_dofs - 1, i, c(i));
    A_aux.AddToEntry(i, N_dofs - 1, c(i));
  }
  /* SAM_LISTING_END_2 */
  /* END_SOLUTION */

  A = A_aux.makeSparse();
  return std::make_pair(A, rhs_vec);
}

}  // namespace RegularizedNeumann

#endif  // define __GRADPROJECTION_H
