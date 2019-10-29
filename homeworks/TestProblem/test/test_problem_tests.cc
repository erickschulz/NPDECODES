#include <gtest/gtest.h>
#include "../mysolution/test_problem.h"

namespace TestProblem {

TEST(TestProblem, TestQuestion1aAdd) {
  EXPECT_EQ(add(2, 3), 5);
  EXPECT_EQ(add(7, 3), 10);
}

}  // namespace TestProblem
