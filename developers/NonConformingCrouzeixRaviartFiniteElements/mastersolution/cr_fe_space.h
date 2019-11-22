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

/** @brief Crouzeix-Raviart finite element space
 *  This class complies with Lehrfempp lf::uscalfe::ScalarUniformFESpace */
/* SAM_LISTING_BEGIN_1 */
class CRFeSpace : public lf::uscalfe::UniformScalarFESpace<scalar_type> {
 public:
  /** @brief no default constructors*/
  CRFeSpace() = delete;
  CRFeSpace(const CRFeSpace &) = delete;
  CRFeSpace(CRFeSpace &&) noexcept = default;
  CRFeSpace &operator=(const CRFeSpace &) = delete;
  CRFeSpace &operator=(CRFeSpace &&) noexcept = default;

  /** @brief main constructor that sets up the local-to-global index mapping
   * @param mesh_p shared pointer to underlying mesh (immutable) */
  explicit CRFeSpace(std::shared_ptr<const lf::mesh::Mesh> mesh_p)
      : lf::uscalfe::UniformScalarFESpace<scalar_type>(
            std::move(mesh_p), std::make_shared<CRReferenceFiniteElement>(),
            nullptr, nullptr) {}

  ~CRFeSpace() override = default;
};
/* SAM_LISTING_END_1 */

}  // namespace NonConformingCrouzeixRaviartFiniteElements

#endif  // NUMPDE_CRFESPACE_H
