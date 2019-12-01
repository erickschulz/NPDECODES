/**
  * @file mesh.h
  * @brief NPDE homework ElementMatrixComputation code
  * @author Janik Schüttler, edited by Oliver Rietmann
  * @date 03.03.2019
  * @copyright Developed at ETH Zurich
  */

#ifndef MESH_H
#define MESH_H

#include <memory>

namespace lf { namespace mesh { class Mesh; } }

std::shared_ptr<lf::mesh::Mesh> Generate2DTestMesh();

#endif // MESH_H
