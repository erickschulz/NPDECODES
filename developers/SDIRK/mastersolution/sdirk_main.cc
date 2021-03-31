#include <iostream>

#include "sdirk.h"

int main() {
  double rate = SDIRK::cvgSDIRK();
  std::cout << std::endl << "The rate is " << rate << std::endl;
  return 0;
}
