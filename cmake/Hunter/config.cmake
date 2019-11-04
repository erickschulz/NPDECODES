# This file specifies additional data for the dependencies that are imported via hunter

hunter_config(Boost VERSION 1.66.0)

hunter_config(Eigen VERSION 3.3.5)

#hunter_config(lehrfempp VERSION 0.7.20)
hunter_config(lehrfempp
  # Currently using commit 'rename No->Num(#157)'
  URL "https://github.com/craffael/lehrfempp/archive/a3042d5ade52ca8b50c1f25d71bc88a283bdb6e8.tar.gz"
  SHA1 "0b1ab3e0ba01ef71817808b2c008dbc404f67a21"
)

hunter_config(GTest VERSION 1.8.0-hunter-p11)
