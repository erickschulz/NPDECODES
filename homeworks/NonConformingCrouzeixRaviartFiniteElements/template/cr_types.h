/**
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss
 * @date   16.03.2019
 * @copyright Developed at ETH Zurich
 */

#ifndef NUMPDE_CRTYPES_H
#define NUMPDE_CRTYPES_H

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>

namespace NonConformingCrouzeixRaviartFiniteElements {

/** Type for indices into global matrices/vectors */
using gdof_idx_t = lf::assemble::gdof_idx_t;
/** Type for global index of entities */
using glb_idx_t = lf::assemble::glb_idx_t;
/** Type for scalar values in Crouzeix-Raviart space */
using scalar_type = double;
/** Type for vector length/matrix size */
using size_type = lf::assemble::size_type;
/** Type for (co-)dimensions */
using dim_t = lf::assemble::dim_t;
/** Type for indexing sub-entities */
using sub_idx_t = lf::base::sub_idx_t;

}  // namespace NonConformingCrouzeixRaviartFiniteElements

#endif  // NUMPDE_CRTYPES_H
