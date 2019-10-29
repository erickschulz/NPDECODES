/**
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss
 * @date   16.03.2019
 * @copyright Developed at ETH Zurich
 */

#ifndef NUMPDE_CRFESPACE_H
#define NUMPDE_CRFESPACE_H

#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>
#include "cr_reference_finite_element.h"
#include "cr_types.h"

namespace NonConformingCrouzeixRaviartFiniteElements {

/**
 * @brief Crouzeix-Raviart Finite Element space
 */
/* SAM_LISTING_BEGIN_1 */
class CRFeSpace : public lf::uscalfe::UniformScalarFESpace<scalar_type> {
 public:
  /**
   * @brief no default constructors
   */
  CRFeSpace() = delete;
  CRFeSpace(const CRFeSpace &) = delete;
  CRFeSpace(CRFeSpace &&) noexcept = default;
  CRFeSpace &operator=(const CRFeSpace &) = delete;
  CRFeSpace &operator=(CRFeSpace &&) noexcept = default;

  /**
   * @brief main constructor that sets up the local-to-global index mapping
   * @param mesh_ptr shared pointer to underlying mesh (immutable)
   */
  explicit CRFeSpace(std::shared_ptr<const lf::mesh::Mesh> mesh_ptr)
      : lf::uscalfe::UniformScalarFESpace<scalar_type>(
            std::move(mesh_ptr), std::make_shared<CRReferenceFiniteElement>(),
            nullptr, nullptr) {
    // TODO: task 2-14.t)
    /* BEGIN_SOLUTION */
    /* TODO Your implementation goes here! */
    /* END_SOLUTION */
  }

  ~CRFeSpace() override = default;
};
/* SAM_LISTING_END_1 */

}  // namespace NonConformingCrouzeixRaviartFiniteElements

#endif  // NUMPDE_CRFESPACE_H
