# This file specifies additional data for the dependencies that are imported via hunter
hunter_config(Eigen VERSION 3.3.5)
hunter_config(Boost VERSION 1.66.0)

# hunter_config(lehrfempp VERSION 0.7.21)
hunter_config(lehrfempp
  # Currently using commit 'Resolve Issue #174'
  URL "https://github.com/craffael/lehrfempp/archive/cc8824bb827d1ca3e3157b4a3b5b50ecf468216b.tar.gz"
  SHA1 "7371d38baa042dfda4da2b112cc5b00b1860c51c"
)

hunter_config(GTest VERSION 1.10.0)
