/**
 * @file contourplot.cc
 * @brief NPDE homework ContourPlot code
 * @author Oliver Rietmann
 * @date 25.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <fstream>
#include <iostream>

#include "contourplot.h"

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

constexpr double Square(double x) { return x * x; }
constexpr double Cube(double x) { return x * x; }

int main() {
  // Compute isolines using gradF and computeIsolinePoints(...)
  Eigen::MatrixXd isolinePoints = ContourPlot::crookedEgg();

  // Compute isolines using F and computeIsolinePointsDQ(...)
  auto F = [] (Eigen::Vector2d x) {
    return Square(x.squaredNorm()) - Cube(x[0]) - Cube(x[1]);
  };
  Eigen::Vector2d y0(1.5, 0.0);
  double T = 10.0;
  Eigen::MatrixXd isolinePointsDQ = ContourPlot::computeIsolinePointsDQ(F, y0, T);

  // Write .csv file with 4 columns containg the x and y coordinates of the two isolines
  std::ofstream csvFile;
  csvFile.open("contourplot.csv");
  csvFile << isolinePoints.format(CSVFormat) << std::endl;
  csvFile << isolinePointsDQ.format(CSVFormat) << std::endl;
  csvFile.close();
  std::cout << "Generated " CURRENT_BINARY_DIR "/contourplot.csv" << std::endl;

  // Call the Python script to plot the isolines
  std::system("python3 " CURRENT_SOURCE_DIR "/contourplot.py " CURRENT_BINARY_DIR "/contourplot.csv " CURRENT_BINARY_DIR "/contourplot.eps");

  return 0;
}
