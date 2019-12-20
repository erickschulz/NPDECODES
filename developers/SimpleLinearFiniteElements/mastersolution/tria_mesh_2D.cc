#include <cassert>
#include <fstream>
#include <iostream>

#include "tria_mesh_2D.h"

namespace SimpleLinearFiniteElements
{

/**
 * @brief read in mesh file into a costum data structure for this exercise
 * @param filename name of the mesh file
 */
TriaMesh2D::TriaMesh2D(std::string filename)
{
  std::ifstream mesh_file(filename, std::ifstream::in);
  std::cout << "Load mesh " << filename << std::endl;
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

Eigen::Matrix<double, 2, 3> TriaMesh2D::operator[] (int i) const {
  Eigen::Matrix<double, 2, 3> triangle;
  for (int k = 0; k < 3; ++k) {
    triangle.col(k) = Coordinates.row(Elements(i,k));
  }
  return triangle;
}

/**
 * @brief Saves a 3D mesh file for plotting function on the 2D mesh.
 * The new z component contains the values of the function on the 2D mesh.
 * @param filename Output file name of the new 3D mesh
 * @param z vector of z values, in correct order
 */
void TriaMesh2D::SaveMesh3D(std::string filename, const Eigen::VectorXd &z) const {

  int n_vertices = Coordinates.rows();
  int n_elements = Elements.rows();

  Eigen::MatrixXd new_vertices(n_vertices, 3);
  new_vertices.leftCols<2>() = Coordinates;
  new_vertices.col(2) = z;

  std::ofstream file(filename);
  if (file.is_open()) {
    file << n_vertices << " Vertices" << std::endl;
    file << new_vertices << std::endl;

    file << n_elements << " Elements" << std::endl;
    file << Elements;
  }
  file.close();
}

} // namespace SimpleLinearFiniteElements
