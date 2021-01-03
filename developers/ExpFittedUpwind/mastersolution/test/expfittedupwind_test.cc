#include <gtest/gtest.h>

#include <string>

#include <Eigen/Core>

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

#include "../expfittedupwind.h"

namespace ExpFittedUpwind::test {

TEST(Bernoulli, Monotinicity) {
    std::vector<double> values = {3.0 , 1.5, 0.5 ,1.0E-5, 1.0E-7, 1.0E-9, -1.0E-9, -1.0E-7, -1.0E-5, -1.0, -5.0, -6.0};

    for(int i = 0; i< values.size()-1; ++i){
        EXPECT_LE(Bernoulli(values[i]), Bernoulli(values[i+1]));
    }
}

TEST(Bernoulli, Cancellation){
    for(int i = 8.0; i >= 10E-12; i/=2.0){
        EXPECT_LE(Bernoulli(i), 2.0);
    }
}


}  // namespace OutputImpedanceBVP::test

