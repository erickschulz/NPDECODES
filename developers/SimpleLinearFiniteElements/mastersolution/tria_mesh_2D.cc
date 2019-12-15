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
  Vertices.reserve(nVertices);
  int nV = 0;
  double tmp1;
  double tmp2;
  while (nV < nVertices)
  {
    mesh_file >> tmp1;
	mesh_file >> tmp2;
	Eigen::Vector2d tmp; 
	tmp << tmp1, tmp2;
	Vertices.push_back(tmp);
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
  Elements.reserve(nElements);
  int nE = 0;
  int tmp3, tmp4, tmp5;
  while (nE < nElements)
  {
    mesh_file >> tmp3;
    mesh_file >> tmp4;
    mesh_file >> tmp5;
	Eigen::Vector3i tmp;
	tmp << tmp3, tmp4, tmp5;
    Elements.push_back(tmp);
	nE++;
  }
  mesh_file.close();
}

/**
 * @brief Adds a z component to the mesh file
 * @param input_file Filename of the mesh to read from
 * @param input_file Filename of the new mesh with z component
 * @param z vector of z-values, in correct order
 */
void TriaMesh2D::addZComponent(std::string input_file, std::string output_file, const Eigen::VectorXd &z) {
  std::ifstream input;
  input.open(input_file);

  std::ofstream output;
  output.open(output_file);

  int n_vertices;
  std::string dummy1;
  input >> n_vertices >> dummy1;
  if (n_vertices != z.size()) {
    std::cout << "Error: Number of vertices of input file and z need to agree!" << std::endl;
    output.close();
    input.close();
    assert(false);
  }

  output << n_vertices << " " << dummy1 << std::endl;
  
  for (int i = 0; i < n_vertices; ++i) {
    double x, y;
    input >> x >> y;
    output << x << " " << y << " " << z(i) << std::endl;
  }

  int n_elements;
  std::string dummy2;
  input >> n_elements >> dummy2;
  output << n_elements << " " << dummy2 << std::endl;

  for (int i = 0; i < n_elements; ++i) {
    int a, b, c;
    input >> a >> b >> c;
    output << a << " " << b << " " << c << std::endl;
  }

  output.close();
  input.close();
}

} // namespace SimpleLinearFiniteElements
