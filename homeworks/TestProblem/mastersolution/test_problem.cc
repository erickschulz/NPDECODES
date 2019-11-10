#include "test_problem.h"

#include <string>

int main() {
  // use implementation of add function defined in utility.cc
  int a = 5;
  int b = 2;
  std::cout << TestProblem::add(a, b);

  // demo reading and writing mesh mesh

  // get path to mesh
  std::string mesh_file = CURRENT_SOURCE_DIR"/meshes/square.msh";

  // read mesh
  std::cout << "Reading mesh from file " << mesh_file << std::endl;
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();

  std::cout << "Shape regularity measure = "
            << TestProblem::shapeRegularityMeasure(*mesh_p) << std::endl;

  // Python output
  lf::io::writeMatplotlib(*mesh_p, "test_mesh.csv");

  // VTK output
  lf::io::VtkWriter vtk_writer(mesh_p, "test_mesh.vtk");

  return 0;
}
