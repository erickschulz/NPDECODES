#include <iostream>
#include <memory>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <lf/mesh/mesh.h>

#include "incidencematrices.h"

int main() {
  // Step 0: Create demo mesh
  std::cout << "Creating demo mesh from exercise sheet.. ";
  std::shared_ptr<lf::mesh::Mesh> demoMesh =
      IncidenceMatrices::createDemoMesh();
  std::cout << "Done!\n\n";

  // Step 1: Compute and print matrix G
  std::cout << "Computing matrix G.. ";
  Eigen::SparseMatrix<int> G =
      IncidenceMatrices::computeEdgeVertexIncidenceMatrix(*demoMesh);
  std::cout << "Done!\n";
  std::cout << "G = \n" << Eigen::MatrixXi(G) << "\n\n";

  // Step 2: Compute and print matrix D
  std::cout << "Computing matrix D.. ";
  Eigen::SparseMatrix<int> D =
      IncidenceMatrices::computeCellEdgeIncidenceMatrix(*demoMesh);
  std::cout << "Done!\n";
  std::cout << "D = \n" << Eigen::MatrixXi(D) << "\n\n";

  // Step 3: Test co-chain complex property (D*G = 0?)
  std::cout << "D*G = 0? "
            << (IncidenceMatrices::testZeroIncidenceMatrixProduct(*demoMesh)
                    ? "Yes!"
                    : "No!")
            << "\n";
  return 0;
}
