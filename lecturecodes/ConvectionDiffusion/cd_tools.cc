/**
 * @file cd_tools.cc
 * @brief Utility Functions for the Convection-Diffusion Problem
 * @author Philippe Peter
 * @date July 2021
 * @copyright Developed at SAM, ETH Zurich
 */

#include "cd_tools.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <algorithm>
#include <memory>

namespace ConvectionDiffusion {

double Diameter(const lf::mesh::Entity& entity) {
  const lf::geometry::Geometry* geo_p = entity.Geometry();
  Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

  switch (entity.RefEl()) {
    case lf::base::RefEl::kTria(): {
      // Diameter of a triangle corresponds to the longest edge
      Eigen::Vector2d e0 = corners.col(1) - corners.col(0);
      Eigen::Vector2d e1 = corners.col(2) - corners.col(1);
      Eigen::Vector2d e2 = corners.col(0) - corners.col(1);
      return std::max(e0.norm(), std::max(e1.norm(), e2.norm()));
    }
    case lf::base::RefEl::kQuad(): {
      // Diameter of a (convex) quadrilateral corresponds to the longer diagonal
      Eigen::Vector2d d0 = corners.col(2) - corners.col(0);
      Eigen::Vector2d d1 = corners.col(3) - corners.col(1);
      return std::max(d0.norm(), d1.norm());
    }
    default: {
      LF_ASSERT_MSG(false,
                    "Diameter not available for " << entity.RefEl().ToString());
    }
  }
}  // namespace ConvectionDiffusion

double MeshWidth(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  double h = 0.0;
  for (const lf::mesh::Entity* entity_p : mesh_p->Entities(0)) {
    h = std::max(h, Diameter(*entity_p));
  }
  return h;
}

}  // namespace ConvectionDiffusion
