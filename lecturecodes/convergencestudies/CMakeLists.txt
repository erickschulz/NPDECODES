add_custom_target(convergencestudies_run)


add_executable(lecturecodes.convergencestudies.twopointbvp twopointbvp.cc)
target_link_libraries(lecturecodes.convergencestudies.twopointbvp
    PUBLIC Eigen3::Eigen
    Boost::program_options
    LF::lf.quad
)
add_custom_target(lecturecodes.convergencestudies.twopointbvp_run
    COMMAND lecturecodes.convergencestudies.twopointbvp -o ${CMAKE_CURRENT_BINARY_DIR}/results_twopointbvp.csv && python3 ${CMAKE_CURRENT_SOURCE_DIR}/twopointbvp_plot.py ${CMAKE_CURRENT_BINARY_DIR}/results_twopointbvp.csv
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/twopointbvp_plot.py
)
add_dependencies(convergencestudies_run lecturecodes.convergencestudies.twopointbvp_run)


add_executable(lecturecodes.convergencestudies.asymptotic asymptotic.cc)
target_link_libraries(lecturecodes.convergencestudies.asymptotic
    PUBLIC Eigen3::Eigen
    Boost::program_options
    LF::lf.quad
)
add_custom_target(lecturecodes.convergencestudies.asymptotic_run
    COMMAND lecturecodes.convergencestudies.asymptotic -o ${CMAKE_CURRENT_BINARY_DIR}/results_asymptotic.csv && python3 ${CMAKE_CURRENT_SOURCE_DIR}/asymptotic_plot.py ${CMAKE_CURRENT_BINARY_DIR}/results_asymptotic.csv
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/asymptotic_plot.py
)
add_dependencies(convergencestudies_run lecturecodes.convergencestudies.asymptotic_run)
