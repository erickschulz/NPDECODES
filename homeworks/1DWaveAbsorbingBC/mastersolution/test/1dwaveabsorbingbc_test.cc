/**
 * @file 1dwaveabsorbingbc_test_mastersolution.cc
 * @brief NPDE homework "1DWaveAbsorbingBC" code
 * @author Oliver Rietmann
 * @date 08.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "../1dwaveabsorbingbc.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace WaveAbsorbingBC1D::test {

TEST(WaveAbsorbingBC1D, waveLeapfrogABC) {
  double c = 1.0;
  double T = 2.0;
  unsigned int N = 4;
  unsigned int m = 10;

  // to test
  Eigen::MatrixXd R = waveLeapfrogABC(c, T, N, m);

  // reference
  Eigen::MatrixXd R_ref(m + 1, N + 1);
  R_ref << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.198669330795061, 0, 0, 0,
      0.127148371708839, 0.389418342308651, 0, 0, 0.0813749578936571,
      0.340774566707901, 0.564642473395035, 0, 0.0520799730519405,
      0.276685692376489, 0.531660472345612, 0.717356090899523,
      0.0370346475036021, 0.21457642371835, 0.491432625671849,
      0.678207714677588, 0.841470984807897, 0.167401538090046,
      0.440634106857561, 0.648527646680643, 0.809707392929288,
      0.932039085967226, 0.376185463723217, 0.624874811472334,
      0.775725839801797, 0.916344317126136, 0.98544972998846, 0.57622943608183,
      0.746498991658528, 0.896375200279672, 0.977212880067294,
      0.999573603041505, 0.719537116976217, 0.855071429793168,
      0.972839902304093, 1.00065619064787, 0.973847630878195, 0.831840148412092,
      0.952273730131951, 0.991735206461539, 0.989139598435836,
      0.909297426825682;

  double error = (R - R_ref).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-8;
  EXPECT_NEAR(error, 0.0, tol);
}

TEST(WaveAbsorbingBC1D, computeEnergies) {
  double c = 1.0;
  double T = 2.0;
  unsigned int N = 4;
  unsigned int m = 10;

  // to test
  Eigen::MatrixXd R = waveLeapfrogABC(c, T, N, m);
  std::pair<Eigen::VectorXd, Eigen::VectorXd> energies =
      computeEnergies(R, c, T / m);
  Eigen::VectorXd E_pot = energies.first;
  Eigen::VectorXd E_kin = energies.second;

  // reference
  Eigen::VectorXd E_pot_ref(m + 1);
  E_pot_ref << 0, 0.0789390059971149, 0.16990449181327, 0.248053760939314,
      0.305310107833622, 0.339420536216696, 0.317639428515545,
      0.218303077367652, 0.116978663744621, 0.0674628156849356,
      0.044885930536866;

  Eigen::VectorXd E_kin_ref(m);
  E_kin_ref << 0.061671098435246, 0.107372816015359, 0.211280571935815,
      0.277989875169054, 0.319954600240159, 0.33022543958348, 0.264741135691184,
      0.166131964617881, 0.0899496009771822, 0.0572728170718303;

  double tol = 1.0e-8;
  double E_pot_error = (E_pot - E_pot_ref).lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(E_pot_error, 0.0, tol);
  double E_kin_error = (E_kin - E_kin_ref).lpNorm<Eigen::Infinity>();
  EXPECT_NEAR(E_kin_error, 0.0, tol);
}

}  // namespace WaveAbsorbingBC1D::test
