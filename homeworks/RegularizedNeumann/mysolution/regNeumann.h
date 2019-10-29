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
    const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> &fe_space,
    const FUNCT_F &f, const FUNCT_H &h) {
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

  // Now fix the solution to be 0 at the node p
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  A = A_aux.makeSparse();
  return std::make_pair(A, rhs_vec);
}

// SUB-EXERCISE f)

// You can use write a helper class which should implement ENTITY_VECTOR_PROVIDER to calculate vector c using
// LehrFEM assembly functions
class VecHelper {
 public:
  explicit VecHelper() {}
  bool isActive(const lf::mesh::Entity &entity) const { return true; }
  Eigen::Vector3d Eval(const lf::mesh::Entity &entity) {

    LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                  "Function only defined for triangular cells");
    Eigen::Vector3d result;
    /* BEGIN_SOLUTION */
    /* TODO Your implementation goes here! */
    /* END_SOLUTION */
    return result;
  }
};

Eigen::VectorXd assembleVector_c( const lf::assemble::DofHandler &dofh){

  Eigen::VectorXd c(dofh.NoDofs());
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return c;

}

template <typename FUNCT_F, typename FUNCT_H>
std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> getGalerkinLSE_augment(
    const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> &fe_space,
    const FUNCT_F &f, const FUNCT_H &h) {
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
  // co-dimension 1 because we locally assemble on edges
  lf::assemble::AssembleVectorLocally(1, dofh, my_vec_provider_edge, rhs_vec);

  // Calculate c and insert it into A
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */

  A = A_aux.makeSparse();
  return std::make_pair(A, rhs_vec);
}

}  // namespace RegularizedNeumann

#endif  // define __GRADPROJECTION_H
