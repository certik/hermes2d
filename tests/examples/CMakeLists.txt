find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})

# examples
add_subdirectory(newton-np-timedep-adapt-system)
add_subdirectory(iron-water)
add_subdirectory(saphir)
add_subdirectory(bracket)
add_subdirectory(bracket-2)
add_subdirectory(singular-perturbation)
