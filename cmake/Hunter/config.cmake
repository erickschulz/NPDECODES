# This file specifies additional data for the dependencies that are imported via hunter
hunter_config(Eigen VERSION 3.3.5)
hunter_config(Boost VERSION 1.66.0)

hunter_config(lehrfempp VERSION 0.7.21)
hunter_config(lehrfempp
  # Currently using commit 'Update README.md'
  URL "https://github.com/craffael/lehrfempp/archive/d36d54238ac0035243a924a3988642b70614d7ba.tar.gz"
  SHA1 "a113e5315e2b1c588684d3e812f5a3adc8954f5a"
)

hunter_config(GTest VERSION 1.10.0)
