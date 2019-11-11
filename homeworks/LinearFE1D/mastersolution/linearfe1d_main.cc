/**
 * @ file LinearFE1D.cc
 * @ brief NPDE homework LinearFE1D code
 * @ author Christian Mitsch, Amélie Loher
 * @ date 11.11.2019
 * @ copyright Developed at ETH Zurich
 */

/* This problem does not rely on Lehrfempp. It is solely based on the
library Eigen. As the computational domain is simply the 1D interval
(0,1), the mesh is stored as a vector whose entries are the
coordinates of the nodes (i.e. the distance of the node from origin 0.0).
*/

#include "linearfe1d.h"

int main() {
  // SOLVING BVP (A), (B) and (C)
  Eigen::VectorXd mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;
  auto alpha = [](double x) { return x; };
  auto f = [](double x) { return x; };
  auto gamma = [](double x) { return x; };
  Eigen::VectorXd uA, uB, uC;
  uA = LinearFE1D::solveA(mesh, gamma, f);
  uB = LinearFE1D::solveB(mesh, alpha, f, 0.1, 0.5);
  uC = LinearFE1D::solveC(mesh, alpha, gamma);

  // PRINTING RESULTS TO.csv FILE
  // Defining CSV output file format
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream ofs;
  // Printing results to file for problem (A)
  std::string filename = "uA.csv";
  ofs.open(filename.c_str());
  if (ofs.is_open()) {
    ofs << uA.format(CSVFormat);
  }
  ofs.close();
  if (ofs.is_open()) {
    std::cout << "File uA.csv was not properly closed." << std::endl;
  }
  // Printing results to file for problem (B)
  filename = "uB.csv";
  ofs.open(filename.c_str());
  if (ofs.is_open()) {
    ofs << uB.format(CSVFormat);
  }
  ofs.close();
  if (ofs.is_open()) {
    std::cout << "File uB.csv was not properly closed." << std::endl;
  }
  // Printing results to file for problem (C)
  filename = "uC.csv";
  ofs.open(filename.c_str());
  if (ofs.is_open()) {
    ofs << uC.format(CSVFormat);
  }
  ofs.close();
}
