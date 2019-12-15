/**
 * @file
 * @brief NPDE homework SimpleLinearFiniteElements
 * @author Amélie Loher
 * @date 11/12/2019
 * @copyright Developed at ETH Zurich
 */

#include <string>
#include <vector>

#include <Eigen/Core>

namespace SimpleLinearFiniteElements {

/**
 * @brief simple mesh data structure used for this problem
 */
struct TriaMesh2D {
  // Constructor: reads mesh data from file
  TriaMesh2D(std::string filename);
  virtual ~TriaMesh2D(void) {}

  static void addZComponent(std::string input_file, std::string output_file, const Eigen::VectorXd &z);

  // Data members describing geometry and topolgy
  std::vector<Eigen::Vector2d> Vertices;
  std::vector<Eigen::Vector3i> Elements;

};

} // namespace SimpleLinearFiniteElements
