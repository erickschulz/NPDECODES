#include <gtest/gtest.h>

#include <string>

#include <Eigen/Core>

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

#include "../expfittedupwind.h"

namespace ExpFittedUpwind::test {

TEST(Bernoulli, Cancellation) {
    std::vector<double> values = {3.0 , 1.5, 0.5 ,1.0E-5, 1.0E-7, 1.0E-9  };
}

}  // namespace OutputImpedanceBVP::test

