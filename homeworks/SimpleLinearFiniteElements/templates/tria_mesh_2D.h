/**
 * @file
 * @brief NPDE homework SimpleLinearFiniteElements
 * @author Amélie Loher
 * @date 11/12/2019
 * @copyright Developed at ETH Zurich
 */

#ifndef TRIAMESH2D_H_
#define TRIAMESH2D_H_

#include <string>

#include <Eigen/Core>

namespace SimpleLinearFiniteElements {

/**
 * @brief simple mesh data structure used for this problem
 */
struct TriaMesh2D {
  // Constructor: reads mesh data from file
  TriaMesh2D(std::string filename);

  Eigen::Matrix<double, 2, 3> operator[] (int i) const;

  void SaveMesh3D(std::string filename, const Eigen::VectorXd &z) const;

  // Data members describing geometry and topolgy
  Eigen::Matrix<double, Eigen::Dynamic, 2> vertices;
  Eigen::Matrix<int, Eigen::Dynamic, 3> elements;
};

}  // namespace SimpleLinearFiniteElements

#endif  // TRIAMESH2D_H_
