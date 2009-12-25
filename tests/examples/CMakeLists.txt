enable_language(Fortran)
find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})
find_package(UMFPACK REQUIRED)
find_package(BLAS REQUIRED)

# examples
add_subdirectory(iron-water)
add_subdirectory(saphir)
add_subdirectory(newton-np-timedep-adapt-system)
