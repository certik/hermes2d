enable_language(Fortran)
find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})
find_package(UMFPACK REQUIRED)
find_package(BLAS REQUIRED)

# mesh tests
add_subdirectory(convert_to_quads)
add_subdirectory(convert_to_triangles)
add_subdirectory(refinements)
add_subdirectory(copy)
add_subdirectory(loader)

