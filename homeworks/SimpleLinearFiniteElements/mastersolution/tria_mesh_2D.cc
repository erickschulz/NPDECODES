#include "simple_linear_finite_elements.h"

#include <iostream>

#include <boost/filesystem.hpp>

namespace SimpleLinearFiniteElements
{

/**
 * @brief read in mesh file into a costum data structure for this exercise
 * @param filename name of the mesh file
 */
TriaMesh2D::TriaMesh2D(const std::string &filename)
{
  boost::filesystem::path here = __FILE__;
  auto filepath =
      (here.parent_path().parent_path() / ("meshes/" + filename)).string();

  std::ifstream mesh_file(filepath, std::ifstream::in);
  std::cout << "Load mesh " << filepath << std::endl;
  if (!mesh_file.good())
  {
    throw std::runtime_error("Cannot open mesh file! File not found");
    return;
  }
  int nVertices;
  mesh_file >> nVertices;
  std::cout << nVertices << " Vertices" << std::endl;
  char keyword[1024];
  mesh_file.getline(keyword, 1024);
  if (!strcmp(keyword, "Vertices"))
  {
    throw std::runtime_error("Keyword 'Vertices' not found. Wrong file format");
    return;
  }
  Coordinates.resize(nVertices, 2);
  int nV = 0;
  while (nV < nVertices)
  {
    mesh_file >> Coordinates(nV, 0);
    mesh_file >> Coordinates(nV, 1);
    nV++;
  }
  int nElements;
  mesh_file >> nElements;
  mesh_file.getline(keyword, 1024);
  if (!strcmp(keyword, "Elements"))
  {
    throw std::runtime_error("Keyword 'Elements' not found. Wrong file format");
    return;
  }
  Elements.resize(nElements, 3);
  int nE = 0;
  while (nE < nElements)
  {
    mesh_file >> Elements(nE, 0);
    mesh_file >> Elements(nE, 1);
    mesh_file >> Elements(nE, 2);
    nE++;
  }
  mesh_file.close();
}

} // namespace SimpleLinearFiniteElements
