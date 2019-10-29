#include "incidence_mat.h"

using namespace Eigen;
using lf::mesh::Mesh;

int main() {
  // Step 0: Create demo mesh
  std::cout << "Creating demo mesh from exercise sheet.. ";
  std::shared_ptr<Mesh> demoMesh = IncidenceMatrices::createDemoMesh();
  std::cout << "Done!\n\n";

  // Step 1: Compute and print matrix G
  std::cout << "Computing matrix G.. ";
  SparseMatrix<int> G =
      IncidenceMatrices::computeEdgeVertexIncidenceMatrix(*demoMesh);
  std::cout << "Done!\n";
  std::cout << "G = \n" << MatrixXi(G) << "\n\n";

  // Step 2: Compute and print matrix D
  std::cout << "Computing matrix D.. ";
  SparseMatrix<int> D =
      IncidenceMatrices::computeCellEdgeIncidenceMatrix(*demoMesh);
  std::cout << "Done!\n";
  std::cout << "D = \n" << MatrixXi(D) << "\n\n";

  // Step 3: Test co-chain complex property (D*G = 0?)
  std::cout << "D*G = 0? "
            << (IncidenceMatrices::testZeroIncidenceMatrixProduct(*demoMesh)
                    ? "Yes!"
                    : "No!")
            << "\n";
  return 0;
}
