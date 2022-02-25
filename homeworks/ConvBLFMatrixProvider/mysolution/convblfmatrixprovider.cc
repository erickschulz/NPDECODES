/**
 * @file
 * @brief demonstration of assembly of Galerkin linear system in LehrFEM++
 * assemble module; meant to provide sample codes for lecture document
 * @author Ralf Hiptmair
 * @date   April 2021
 * @copyright Developed at ETH Zurich
 */

#include "convblfmatrixprovider.h"

namespace cblfdemo {

std::pair<Eigen::Matrix<double, 2, 3>, double> getGradBaryCoords(
    const lf::mesh::Entity& tria) {
  // Throw error in case no triangular cell
  LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(),
                "getGradBeyCoords: Unsupported cell type " << tria.RefEl());
  // Obtain vertex coordinates of the triangle in a 2x3 matrix
  const auto vertices = lf::geometry::Corners(*(tria.Geometry()));
  // Set up an auxiliary 3x3-matrix with a leading column 1 and
  // the vertex coordinates in its right 3x2 block
  Eigen::Matrix<double, 3, 3> X;  // temporary matrix
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  // The determinant of the auxiliary matrix also supplies the determinant
  // Compute the gradients of the barycentric coordinate functions
  // and store them in the columns of a 2x3 matrix
  return {X.inverse().block<2, 3>(1, 0), 0.5 * std::abs(X.determinant())};
}

double testCDBLF(std::shared_ptr<lf::mesh::Mesh> mesh_p, Eigen::Vector2d a,
                 Eigen::Vector2d b) {
  // Build linear mesh functions
  auto vec_a = [a](Eigen::Vector2d /*x*/) -> Eigen::Vector2d { return a; };
  auto lin_b = [b](Eigen::Vector2d x) -> double { return b.dot(x); };
  const lf::mesh::utils::MeshFunctionGlobal MF_va(vec_a);
  const lf::mesh::utils::MeshFunctionGlobal MF_b(lin_b);
  // DofHandler for linear Lagrangian FE space
  lf::assemble::UniformFEDofHandler dh_linfe(mesh_p,
                                             {{lf::base::RefEl::kPoint(), 1}});
  // DofHandler for piecewise constants
  lf::assemble::UniformFEDofHandler dh_constfe(mesh_p,
                                               {{lf::base::RefEl::kTria(), 1}});
  // Matrix in triplet format holding temporary Galerkin matrix
  const unsigned int N_linfe = dh_linfe.NumDofs();
  const unsigned int N_constfe = dh_constfe.NumDofs();
  lf::assemble::COOMatrix<double> mat(N_constfe, N_linfe);
  // Initialize Galerkin matrix
  CDBLFElemMatProvider cdblfelmatprovider(MF_va);
  mat = lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, dh_linfe, dh_constfe, cdblfelmatprovider);

  // Compute coefficent vector of nodal interpolant
  lf::uscalfe::FeSpaceLagrangeO1<double> fes_lin(mesh_p);
  auto cv_interp = lf::fe::NodalProjection(fes_lin, MF_b);

  // Indirectly compute integral of inner product of the two vectors over domain
  const Eigen::SparseMatrix<double> A(mat.makeSparse());
  return Eigen::VectorXd::Constant(N_constfe, 1.0).transpose() * A * cv_interp;
}

}  // namespace cblfdemo
