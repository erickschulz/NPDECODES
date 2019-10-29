#include "test_problem.h"

int main() {
  // use implementation of add function defined in utility.cc
  int a = 5;
  int b = 2;
  std::cout << TestProblem::add(a, b);

  // demo reading and writing mesh mesh

  // get path to mesh
  boost::filesystem::path here = __FILE__;
  auto mesh_path = here.parent_path().parent_path() / "meshes/square.msh";

  // read mesh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_path.string());
  auto mesh = reader.mesh();

  // Python output
  lf::io::writeMatplotlib(*mesh, "test_mesh.csv");

  // VTK output
  lf::io::VtkWriter vtk_writer(mesh, "test_mesh.vtk");

  return 0;
}
