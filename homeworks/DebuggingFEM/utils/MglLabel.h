#ifndef MGL_LABEL_HPP
#define MGL_LABEL_HPP

#include <mgl2/mgl.h>

namespace mgl {

struct MglLabel {
  MglLabel(const std::string &str ="", const double pos=0)
    : str_(str)
    , pos_(pos)
  {}

  std::string str_;
  double pos_;
};

}
#endif
