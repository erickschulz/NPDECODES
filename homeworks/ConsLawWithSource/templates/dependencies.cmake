# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/conslawwithsource_main.cc
  ${DIR}/conslawwithsource.h
  ${DIR}/conslawwithsource.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
