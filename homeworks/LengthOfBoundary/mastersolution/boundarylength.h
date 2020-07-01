#ifndef __BOUNDARYLENGTH_H
#define __BOUNDARYLENGTH_H

/**
 * @ file boundarylength.h
 * @ brief NPDE homework LengthOfBoundary code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include <memory>
#include <string>

#include <lf/mesh/mesh.h>

namespace LengthOfBoundary {

double volumeOfDomain(const std::shared_ptr<lf::mesh::Mesh> mesh_p);

double lengthOfBoundary(const std::shared_ptr<lf::mesh::Mesh> mesh_p);

std::pair<double, double> measureDomain(std::string filename);

} // namespace LengthOfBoundary
#endif
