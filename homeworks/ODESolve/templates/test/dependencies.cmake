# Add your custom dependencies here:

# PROBLEM_NAME and DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/odesolve_test.cc
  ${DIR}/odesolve.h
  ${DIR}/odesolve.cc
  ${DIR}/polyfit.h
  ${DIR}/polyfit.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
)
