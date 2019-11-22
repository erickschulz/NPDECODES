# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/laxwendroffscheme_main.cc
  ${DIR}/laxwendroffscheme.h
  ${DIR}/laxwendroffscheme.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
