#ifndef SYMPLECTIC_ASSEMBLE_HPP
#define SYMPLECTIC_ASSEMBLE_HPP
/**
 * @file
 * @brief NPDE homework SymplecticTimesteppingWaves
 * @author Erick Schulz
 * @date 25/03/2019
 * @copyright Developed at ETH Zurich
 */

// common includes
#include <math.h>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <string>
#include <chrono>
#include <iostream>
#include <thread>

// Eigen and boost includes
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <boost/assert.hpp>
#include <boost/filesystem.hpp>
#include <unsupported/Eigen/KroneckerProduct>

// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

namespace SymplecticTimesteppingWaves {

template <typename FUNC_ALPHA, typename FUNC_GAMMA, typename FUNC_BETA>
Eigen::SparseMatrix<double> assembleGalerkinMatrix(
    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_p,
    FUNC_ALPHA alpha, FUNC_GAMMA gamma, FUNC_BETA beta) {
  Eigen::SparseMatrix<double> galMat;

  /* SOLUTION_BEGIN */
  /* TODO Your implementation goes here! */
  /* SOLUTION_END */
  return galMat;
}

}  // namespace SymplecticTimesteppingWaves

#endif
