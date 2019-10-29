/** @file
This homework "TransformationOfGalerkinMatrices" consists of matrices
transformation (change of Galerkin bases) using the triplet format.
 */

#include <Eigen/Sparse>
#include <vector>

namespace TransformationOfGalerkinMatrices {

// Eigen's built-in type for COO matrix format
typedef Eigen::Triplet<double> triplet_t;

/** @brief Transformation of Galerking matrix in COO format using the change of
 * basis described in the problem statement.
 * @param A "Old" Galerkin matrix in COO format
 * @return triplets describing "New" Galerkin matrix
 */
std::vector<triplet_t> transformCOOmatrix(const std::vector<triplet_t> &A);

}  // namespace TransformationOfGalerkinMatrices
