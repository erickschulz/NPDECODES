/**
 * @file
 * @brief NPDE homework "Handling degrees of freedom (DOFs) in LehrFEM++"
 * @author Julien Gacon
 * @date March 1st, 2019
 * @copyright Developed at ETH Zurich
 */

#include "handling_dofs.h"

namespace HandlingDOFs {

std::array<std::size_t, 3> countEntityDofs(
    const lf::assemble::DofHandler& dofhandler) {
  std::array<std::size_t, 3> entityDofs;
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return entityDofs;
}

std::size_t countBoundaryDofs(const lf::assemble::DofHandler& dofhandler) {
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  // given an entity, bd_flags(entity) = true, if the entity is on the boundary
  lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags(
      lf::mesh::utils::flagEntitiesOnBoundary(mesh));
  std::size_t no_dofs_on_bd = 0;
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return no_dofs_on_bd;
}

double integrateLinearFEFunction(const lf::assemble::DofHandler& dofhandler,
                                 const Eigen::VectorXd& mu) {
  double I = 0;
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return I;
}

double integrateQuadraticFEFunction(const lf::assemble::DofHandler& dofhandler,
                                    const Eigen::VectorXd& mu) {
  double I = 0;
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return I;
}

Eigen::VectorXd convertDOFsLinearQuadratic(
    const lf::assemble::DofHandler& dofh_Linear_FE,
    const lf::assemble::DofHandler& dofh_Quadratic_FE,
    const Eigen::VectorXd& mu) {
  if (dofh_Linear_FE.Mesh() != dofh_Quadratic_FE.Mesh()) {
    throw "Underlying meshes must be the same for both DOF handlers!";
  }
  std::shared_ptr<const lf::mesh::Mesh> mesh =
      dofh_Linear_FE.Mesh();                         // get the mesh
  Eigen::VectorXd zeta(dofh_Quadratic_FE.NoDofs());  // initialise empty zeta
  // safety guard: always set zero if you're not sure to set every entry later
  // on for us this shouldn't be a problem, but just to be sure
  zeta.setZero();

  for (const auto& cell : mesh->Entities(0)) {
    // check if the spaces are actually linear and quadratic
    /* BEGIN_SOLUTION */
    /* TODO Your implementation goes here! */
    /* END_SOLUTION */

    // get the global dof indices of the linear and quadratic FE spaces, note
    // that the vectors obey the LehrFEM++ numbering, which we will make use of
    // lin_dofs will have size 3 for the 3 dofs on the nodes and
    // quad_dofs will have size 6, the first 3 entries being the nodes and
    // the last 3 the edges
    /* BEGIN_SOLUTION */
    /* TODO Your implementation goes here! */
    /* END_SOLUTION */

    // assign the coefficients of mu to the correct entries of zeta, use
    // the previous subproblem 2-9.a
    /* BEGIN_SOLUTION */
    /* TODO Your implementation goes here! */
    /* END_SOLUTION */
  }
  return zeta;
}

}  // namespace HandlingDOFs
