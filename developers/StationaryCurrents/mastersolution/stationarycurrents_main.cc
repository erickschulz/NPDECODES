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
      // Loop through "bentwire0" - "bentwire6"
      std::vector<std::pair<double, double>> fluxes{};
      for (unsigned int l = 0; l <= 6; ++l) {
        basename = "bentwire" + std::to_string(l);
        fluxes.push_back(dmxbc::computePotential(basename));
      }
      double ref_flux = fluxes[6].second;
      for (unsigned int l = 0; l <= 6; ++l) {
        std::cout << "Level " << l << ": bd flux = " << fluxes[l].first
                  << " (err = " << std::abs(fluxes[l].first - ref_flux)
                  << "), vol flux = " << fluxes[l].second
                  << " (err = " << std::abs(fluxes[l].second - ref_flux) << ")"
                  << std::endl;
      }
      break;
    }
    case 2: {
      // Filename given via the command line
      basename = argv[1];
      auto [bd_flux, vol_flux] = dmxbc::computePotential(basename);
      std::cout << basename << ": fluxes = " << bd_flux << ", " << vol_flux
                << std::endl;
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
