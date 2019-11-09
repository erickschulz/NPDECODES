#include "test_problem.h"

int main() {
  // use implementation of add function defined in utility.cc
  int a = 5;
  int b = 2;
  std::cout << add(a, b);

  // demo reading and writing mesh mesh

  // read mesh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), CURRENT_SOURCE_DIR"/meshes/square.msh");
  auto mesh = reader.mesh();

  // Python output
  lf::io::writeMatplotlib(*mesh, "test_mesh.csv");

  // VTK output
  lf::io::VtkWriter vtk_writer(mesh, "test_mesh.vtk");

  return 0;
}
