#include <lf/uscalfe/uscalfe.h>
#include <Eigen/Dense>

namespace TranspSemiLagr{

template<typename FUNCTOR, typename MESH_FUNCTION>
class UpwindLagrangianElementVectorProvider{
  static_assert(lf::mesh::utils::isMeshFunction<MESH_FUNCTION>);

 public:
  UpwindLagrangianElementVectorProvider(FUNCTOR v,double tau, 
                                        std::shared_ptr<const lf::mesh::Mesh> mesh_p,
                                        MESH_FUNCTION U0 ): v_(v), tau_(tau), mesh_p_(std::move(mesh_p)), U0_(U0) {}
  Eigen::Vector3d Eval(const lf::mesh::Entity& entity);
  bool isActive(const lf::mesh::Entity & /* entity */) const {return true;}

private:

  Eigen::Vector2d pull_back_(Eigen::Vector2d x, const lf::geometry::Geometry& geo);

  double tau_;
  FUNCTOR v_;
  MESH_FUNCTION U0_;
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;
};


template<typename FUNCTOR, typename MESH_FUNCTION>
Eigen::Vector3d
UpwindLagrangianElementVectorProvider<FUNCTOR, MESH_FUNCTION>::Eval(const lf::mesh::Entity& entity){
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  const lf::geometry::Geometry *geo_ptr = entity.Geometry();
  Eigen::Vector3d result;
  result.setZero();

  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
  const double area = lf::geometry::Volume(*geo_ptr);

  for(unsigned int i = 0; i < 3; ++i){
    Eigen::Vector2d y = corners.col(i) - tau_*v_(corners.col(i));

    //iterate over all mesh triangles to find the one containing the evaluation point.
    for( const lf::mesh::Entity* cell: mesh_p_->Entities(0)){
       Eigen::Vector2d yhat = pull_back_(y, *(cell->Geometry()));
        if(yhat(0) >= 0 && yhat(1) >= 0 && yhat.sum() <= 1){
          result(i) = U0_(*cell, yhat)[0] * area / 3.0;
      }
    }
  }
  return result;
}

template<typename FUNCTOR, typename MESH_FUNCTION>
Eigen::Vector2d UpwindLagrangianElementVectorProvider<FUNCTOR, MESH_FUNCTION>::
pull_back_(Eigen::Vector2d x, const lf::geometry::Geometry& geo){
  const Eigen::MatrixXd corners = lf::geometry::Corners(geo);
  const Eigen::MatrixXd InvJacobian = geo.Jacobian(corners.col(0)).inverse();
  return InvJacobian*(x - corners.col(0));
}


template<typename FUNCTOR>
class MassLumpedElementMatrixProvider{
 public:
  MassLumpedElementMatrixProvider(FUNCTOR c): c_(c) {}
  Eigen::Matrix3d Eval(const lf::mesh::Entity& entity);
  bool isActive(const lf::mesh::Entity & /* entity */) const {return true;}

private:
  FUNCTOR c_;

};


template<typename FUNCTOR>
Eigen::Matrix3d MassLumpedElementMatrixProvider<FUNCTOR>::Eval(const lf::mesh::Entity& entity){
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");
  const lf::geometry::Geometry *geo_ptr = entity.Geometry();
  Eigen::Matrix3d result = Eigen::Matrix3d::Zero();


  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
  const double area = lf::geometry::Volume(*geo_ptr);

  for(int i = 0; i < 3; ++i){
    result(i,i) = area / 3.0 * c_(corners.col(i));
  }

  return result;
}

} //namespace TranspSemiLagr
