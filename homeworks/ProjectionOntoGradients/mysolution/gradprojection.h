#ifndef __GRADPROJECTION_H
#define __GRADPROJECTION_H

#include <iostream>
#include <memory>

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace ProjectionOntoGradients {

// Start of sub-problem e)
class ElementMatrixProvider {
 private:
  using coord_t = Eigen::Vector2d;

 public:
  Eigen::Matrix3d Eval(const lf::mesh::Entity &entity);
  bool isActive(const lf::mesh::Entity &entity) const { return true; }
};

Eigen::Matrix3d ElementMatrixProvider::Eval(const lf::mesh::Entity &entity) {
  const lf::geometry::Geometry *geo_ptr = entity.Geometry();
  Eigen::Matrix3d loc_mat;

  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION

  return loc_mat;
}
// End of sub-problem e)

// Start of sub-problem g)
template <typename FUNCTOR>
class GradProjRhsProvider {
 private:
  FUNCTOR f_;
  using coord_t = Eigen::Vector2d;

 public:
  explicit GradProjRhsProvider(FUNCTOR f) : f_(f) {}
  Eigen::Vector3d Eval(const lf::mesh::Entity &entity);
  bool isActive(const lf::mesh::Entity &entity) const { return true; }
};

template <typename FUNCTOR>
Eigen::Vector3d GradProjRhsProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity &entity) {
  Eigen::Vector3d loc_vec;

  const lf::geometry::Geometry *geo_ptr = entity.Geometry();

  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION

  return loc_vec;
}
// End of sub-problem e)

// Start of sub-problem h)
template <typename FUNCTOR>
Eigen::VectorXd projectOntoGradients(const lf::assemble::DofHandler &dofh,
                                     FUNCTOR f) {
  const std::size_t N_dofs = dofh.NoDofs();
  Eigen::VectorXd sol_vec;

  // ASSEMBLE GLOBAL MATRIX
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION
  // ASSEMBLE GLOBAL RHS VECOR
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION
  // ENFORCE HOMOGENEOUS DIRICHLET BOUNDARY CONDITIONS
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION
  // CALCULATE THE SOLUTION
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION
  return sol_vec;
}
// End of sub-problem h)

}  // namespace ProjectionOntoGradients 

#endif  // define __GRADPROJECTION_H
