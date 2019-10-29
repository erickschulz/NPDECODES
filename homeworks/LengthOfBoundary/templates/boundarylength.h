#ifndef __BOUNDARYLENGTH_H
#define __BOUNDARYLENGTH_H

/**
 * @ file boundarylength.h
 * @ brief NPDE homework LengthOfBoundary code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include <iostream>
#include <memory>
#include <string>

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>

namespace LengthOfBoundary {

double volumeOfDomain(const std::shared_ptr<lf::mesh::Mesh> mesh);

double lengthOfBoundary(const std::shared_ptr<lf::mesh::Mesh> mesh);

std::pair<double, double> measureDomain(std::string msh_file_name);

}  // namespace LengthOfBoundary
#endif
