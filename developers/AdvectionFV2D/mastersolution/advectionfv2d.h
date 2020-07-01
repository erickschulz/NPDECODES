#ifndef __ADCECTIONFV2D_H
#define __ADCECTIONFV2D_H
/**
 * @file advectionfv2d.h
 * @brief NPDE homework AdvectionFV2D code
 * @author Philipp Egg
 * @date 21.06.2020
 * @copyright Developed at ETH Zurich
 */

#include <array>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

namespace AdvectionFV2D {

/**
 * @brief Short description of this function.
 *
 * @param x Describes the argument x.
 * @param n Describes the argument n.
 * @return Describes the return value.
 */

Eigen::VectorXd dummyFunction(double x, int n);

std::shared_ptr<
    lf::mesh::utils::CodimMeshDataSet<Eigen::Matrix<double, 2, Eigen::Dynamic>>>
computeCellNormals(std::shared_ptr<const lf::mesh::Mesh> mesh_p);

////////////////////////////////////////////////////////////////////////
// TODO: Inconsistenent file name in solution (getAdjacentDofIndex())
////////////////////////////////////////////////////////////////////////
std::shared_ptr<
    lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity *, 4>>>
getAdjacentCellPointers(std::shared_ptr<const lf::mesh::Mesh> mesh_p);

Eigen::Matrix<double, 2, 3> gradbarycoordinates(
    const Eigen::Matrix<double, 2, 3> &triangle);

Eigen::Vector2d barycenter(const Eigen::MatrixXd corners);

template <typename VECTORFIELD>
Eigen::SparseMatrix<double> initializeMOLODEMatrix(
    const lf::assemble::DofHandler &dofh, VECTORFIELD &&beta) {
  // Compute cell normals
  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
      Eigen::Matrix<double, 2, Eigen::Dynamic>>>
      normal_vectors = AdvectionFV2D::computeCellNormals(dofh.Mesh());

  // Compute adjecent cells
  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
      std::array<const lf::mesh::Entity *, 4>>>
      adjacentCells = AdvectionFV2D::getAdjacentCellPointers(dofh.Mesh());

  // Set up matrix B
  int bound_nnz = dofh.Mesh()->NumEntities(0) + dofh.Mesh()->NumEntities(1);
  int num_dof = dofh.NumDofs();
  Eigen::SparseMatrix<double> B_Matrix(num_dof, num_dof);
  B_Matrix.reserve(bound_nnz);

  for (const lf::mesh::Entity *cell : dofh.Mesh()->Entities(0)) {
    // Compute area of cell
    const lf::geometry::Geometry *geo_p = cell->Geometry();
    double area = lf::geometry::Volume(*geo_p);

    // Corresponting DOF of ref. cell
    int row = dofh.GlobalDofIndices(*cell)[0];

    // Corresponting normal vectors of ref. cell
    Eigen::Matrix<double, 2, Eigen::Dynamic> cur_normal_vectors =
        (*normal_vectors)(*cell);

    int num_edges = cur_normal_vectors.cols();
    // Get adjacent cells of ref. cell
    std::array<const lf::mesh::Entity *, 4> neighbour_cells =
        (*adjacentCells)(*cell);

    // Get edges of ref. cell
    auto cell_edges = cell->SubEntities(1);

    int counter = 0;
    for (const lf::mesh::Entity *next_cell : neighbour_cells) {
      // Check whether cell has a neighbor
      if (next_cell != nullptr) {
        // Geo pointer of edge to next_cell
        const lf::geometry::Geometry *geo_p = (cell_edges[counter])->Geometry();

        // Length of current edge
        double edge_length = lf::geometry::Volume(*geo_p);

        // Corners of current edge
        Eigen::Matrix<double, 2, Eigen::Dynamic> corner_edges =
            lf::geometry::Corners(*geo_p);

        // Midpoint of current edge
        Eigen::Vector2d midpoint =
            (corner_edges.col(0) + corner_edges.col(1)) * 0.5;

        // Corresponting DOF of next_cell
        int col = dofh.GlobalDofIndices(*next_cell)[0];

        // Compute Flux and store it directly in B
        double flux = (cur_normal_vectors.col(counter)).dot(beta(midpoint));
        if (flux >= 0) {
          B_Matrix.coeffRef(row, row) -= flux * edge_length / area;
        } else {  // flux < 0
          B_Matrix.coeffRef(row, col) -= flux * edge_length / area;
        }
      } else if (counter < num_edges) {
        // In the case a cell has no neighbor, the flux of the edge has
        // to be considered if it is greater than 0

        // Geo pointer of edge to next_cell
        const lf::geometry::Geometry *geo_p = (cell_edges[counter])->Geometry();

        // Length of current edge
        double edge_length = lf::geometry::Volume(*geo_p);

        // Corners of current edge
        Eigen::Matrix<double, 2, Eigen::Dynamic> corner_edges =
            lf::geometry::Corners(*geo_p);

        // Midpoint of current edge
        Eigen::Vector2d midpoint =
            (corner_edges.col(0) + corner_edges.col(1)) * 0.5;

        // Compute Flux and store it directly in B
        double flux = cur_normal_vectors.col(counter).dot(beta(midpoint));
        if (flux >= 0) {
          B_Matrix.coeffRef(row, row) -= flux * edge_length / area;
        }
      }
      counter++;
    }
  }
  return B_Matrix;
}

////////////////////////////////////////////////////////////////////////
// TODO: Inconsistency - Typo in computeMin arguments (constd -> const)
////////////////////////////////////////////////////////////////////////
double computeHmin(std::shared_ptr<const lf::mesh::Mesh> mesh_p);

template <typename VECTORFIELD>
Eigen::VectorXd solveAdvection2D(const lf::assemble::DofHandler &dofh,
                                 VECTORFIELD &&beta, const Eigen::VectorXd u0_h,
                                 double T, unsigned int M) {
  double tau = T / M;

  // Initialize vector for the result
  int num_dof = dofh.NumDofs();
  Eigen::VectorXd result(num_dof);

  // Compute B
  Eigen::MatrixXd B_matrix = initializeMOLODEMatrix(dofh, beta);

  // Set mu to inital bump
  Eigen::VectorXd mu = u0_h;

  // Enforce Dirichlet Boundary Condition
  // Might interfere with the initial bump!
  if (true) {
    // Flag Edges on the boundary
    lf::mesh::utils::CodimMeshDataSet<bool> bd_blags{
        lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 1)};

    // Get cell normals
    std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<
        Eigen::Matrix<double, 2, Eigen::Dynamic>>>
        normal_vectors = AdvectionFV2D::computeCellNormals(dofh.Mesh());

    // Iterate over all cells
    for (const lf::mesh::Entity *cell : dofh.Mesh()->Entities(0)) {
      auto edges = cell->SubEntities(1);
      int counter = 0;
      // Iterate over all edges
      for (auto edge : edges) {
        // Cell contains edge at boundary
        if (bd_blags(*edge)) {
          const lf::geometry::Geometry *geo_p = edge->Geometry();

          // Corners of current edge
          Eigen::Matrix<double, 2, Eigen::Dynamic> corner_edges =
              lf::geometry::Corners(*geo_p);

          // Midpoint of current edge
          Eigen::Vector2d midpoint =
              (corner_edges.col(0) + corner_edges.col(1)) * 0.5;

          // Corresponting normal vectors of ref. cell
          Eigen::Matrix<double, 2, Eigen::Dynamic> cur_normal_vectors =
              (*normal_vectors)(*cell);

          // If the flux is < 0, set entity in mu and the row of M to zero
          if ((cur_normal_vectors.col(counter)).dot(beta(midpoint)) < 0) {
            int idx = dofh.GlobalDofIndices(*cell)[0];
            B_matrix.row(idx).setZero();
            mu[idx] = 0.0;
          }
        }
      }
      counter++;
    }
  }

  Eigen::VectorXd k0;
  Eigen::VectorXd k1;

  // Timestepping
  for (int step = 0; step < M; step++) {
    k0 = B_matrix * mu;
    k1 = B_matrix * (mu + tau * k0);
    mu = mu + tau * 0.5 * (k0 + k1);

    if (mu.lpNorm<Eigen::Infinity>() > 1000 * u0_h.lpNorm<Eigen::Infinity>()) {
      throw std::overflow_error("Overflow occured!!\n");
    }
  }
  return mu;
}

Eigen::VectorXd simulateAdvection(const lf::assemble::DofHandler &dofh);

Eigen::VectorXd refSolution(const lf::assemble::DofHandler &dofh);

////////////////////////////////////////////////////////////////////////
// TODO: Return int instead of double - OK??
////////////////////////////////////////////////////////////////////////
template <typename VECTORFIELD>
int findCFLthreshold(const lf::assemble::DofHandler &dofh, VECTORFIELD &&beta,
                     double T) {
  // Set upper and lower limit
  int M_upper = int((T / computeHmin(dofh.Mesh())) + 1);
  int M_lower = 1;

  // Initialize a vector for the result and
  // randomly initialize a vector for the initial condition
  int num_dof = dofh.NumDofs();
  Eigen::VectorXd result(num_dof);
  Eigen::VectorXd u0_h = Eigen::VectorXd::Random(num_dof);

  // Shift upper and lower bound until contition below isn't satisfied
  while ((M_upper - M_lower) > 2) {
    // Perform a simulation at M_middle
    // If it succeeds, set the upper bound to M_middle
    // Otherwise set the lower bound to M_middle
    // Repeat ...
    int M_middle = (M_upper + M_lower) / 2;
    try {
      Eigen::VectorXd result = solveAdvection2D(dofh, beta, u0_h, T, M_middle);
      M_upper = M_middle;
    } catch (const std::exception &e) {
      M_lower = M_middle;
    }
  }

  int thres = M_upper;
  std::cout << "Final threshold: " << thres << " vs. "
            << int((1.0 / computeHmin(dofh.Mesh())) + 1) << " (T / hmin)"
            << std::endl;
  return thres;
}

}  // namespace AdvectionFV2D

#endif  // define __ADCECTIONFV2D_H
