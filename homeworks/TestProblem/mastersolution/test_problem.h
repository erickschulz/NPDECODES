#include <iostream>

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <boost/filesystem.hpp>

namespace TestProblem {

/** @brief trivial function */
int add(int a, int b);

/**
 * @brief Computation of shape regularity measure
 *
 * @param mesh a mesh object
 * @return shape regularity measure of mesh
 */
double shapeRegularityMeasure(const lf::mesh::Mesh &mesh);

}  // namespace TestProblem
