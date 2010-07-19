find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})

# examples
add_subdirectory(newton-np-timedep-adapt-system)
add_subdirectory(iron-water)
add_subdirectory(saphir)
add_subdirectory(bracket)
add_subdirectory(singular-perturbation)
add_subdirectory(neutronics-4-group-adapt)
