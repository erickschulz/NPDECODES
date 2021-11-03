/**
 * @ file
 * @ brief NPDE homework problem: Computation of stationary currents
 * @ author Ralf Hiptmair
 * @ date July 2020
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "stationarycurrents.h"

int main(int argc, char** argv) {
  // Expects at most one command line argument, which should be the name of the
  // Gmsh mesh file WITHOUT the .msh suffix.
  std::string basename{};
  switch (argc) {
    case 1: {
      // Loop through meshes "bentwire0" - "bentwire6" for studying convergence
      std::vector<std::tuple<double, double, double>> fluxes{};
      for (unsigned int l = 0; l <= 6; ++l) {
        basename = "bentwire" + std::to_string(l);
        fluxes.push_back(dmxbc::computePotential(basename));
      }
      double ref_flux = 0.1859836202175363;  // "Reference value"
      const int fieldwidth = 10;
      std::cout << std::fixed << std::setprecision(6) << std::setfill(' ');
      std::cout << std::setw(fieldwidth) << "h" << std::setw(fieldwidth)
                << "bd_flux" << std::setw(fieldwidth) << "err(bd)"
                << std::setw(fieldwidth) << "vol_flux" << std::setw(fieldwidth)
                << "err(vol)" << std::endl;
      for (unsigned int l = 0; l <= 6; ++l) {
        const double h = std::get<0>(fluxes[l]);
        const double bd_flux = std::get<1>(fluxes[l]);
        const double vol_flux = std::get<2>(fluxes[l]);
        std::cout << std::setw(fieldwidth) << h << std::setw(fieldwidth)
                  << bd_flux << std::setw(fieldwidth)
                  << std::abs(bd_flux - ref_flux) << std::setw(fieldwidth)
                  << vol_flux << std::setw(fieldwidth)
                  << std::abs(vol_flux - ref_flux) << std::endl;
      }
      break;
    }
    case 2: {
      // Filename given via the command line
      basename = argv[1];
      auto [h, bd_flux, vol_flux] = dmxbc::computePotential(basename);
      std::cout << basename << ": h = " << h << ", fluxes = " << bd_flux << ", "
                << vol_flux << std::endl;
      break;
    }
    default: {
      std::cerr << "Usage: " << argv[0] << " <basename (no suffix!)> "
                << std::endl;
      break;
    }
  }  // end switch

  return 0;
}
