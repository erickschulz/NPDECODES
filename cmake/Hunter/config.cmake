# This file specifies additional data for the dependencies that are imported via hunter

hunter_config(Boost VERSION 1.66.0)

hunter_config(Eigen VERSION 3.3.5)

#hunter_config(lehrfempp VERSION 0.7.20)
hunter_config(lehrfempp
  # Currently using commit 'rename No->Num(#157)'
  URL "https://github.com/craffael/lehrfempp/archive/489d6da16187b44026c2c94f2592daa0dfe7170f.tar.gz"
  SHA1 "3bd5fe8b9013f6edbb784d40fcc8c7b402fca22c"
)

hunter_config(GTest VERSION 1.10.0)
