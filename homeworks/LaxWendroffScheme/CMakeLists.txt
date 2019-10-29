include(GoogleTest)

set(template_sources
  templates/laxwendroffscheme_main.cc
  templates/laxwendroffscheme.h
  templates/laxwendroffscheme.cc
)
set(mastersolution_sources
  mastersolution/laxwendroffscheme_main.cc
  mastersolution/laxwendroffscheme.h
  mastersolution/laxwendroffscheme.cc
)
configure_file(mastersolution/plot.py plot_mastersolution.py COPYONLY)

set(mysolution_sources
  mysolution/laxwendroffscheme_main.cc
  mysolution/laxwendroffscheme.h
  mysolution/laxwendroffscheme.cc
)

set(mastersolution_test_sources
  test/laxwendroffscheme_test_mastersolution.cc
)

set(mysolution_test_sources
  test/laxwendroffscheme_test_mysolution.cc
)

# Libraries
add_library(LaxWendroffScheme_mastersolution.solution ${mastersolution_sources})

add_library(LaxWendroffScheme_mysolution.solution ${mysolution_sources})

# Executables
add_executable(LaxWendroffScheme_template ${template_sources})

add_executable(LaxWendroffScheme_mysolution ${mysolution_sources})

add_executable(LaxWendroffScheme_mastersolution ${mastersolution_sources})

add_executable(LaxWendroffScheme_test_mastersolution ${mastersolution_test_sources})

add_executable(LaxWendroffScheme_test_mysolution ${mysolution_test_sources})


target_link_libraries(LaxWendroffScheme_mastersolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(LaxWendroffScheme_mysolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(LaxWendroffScheme_template
	PUBLIC Eigen3::Eigen
)

target_link_libraries(LaxWendroffScheme_mastersolution.solution
  PUBLIC Eigen3::Eigen 
)

target_link_libraries(LaxWendroffScheme_mysolution.solution
  PUBLIC Eigen3::Eigen 
)

target_link_libraries(LaxWendroffScheme_test_mastersolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	LaxWendroffScheme_mastersolution.solution
)

target_link_libraries(LaxWendroffScheme_test_mysolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	LaxWendroffScheme_mysolution.solution
)

gtest_discover_tests(LaxWendroffScheme_test_mysolution)
gtest_discover_tests(LaxWendroffScheme_test_mastersolution)
