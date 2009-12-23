project(convert_to_triangles)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(convert_to_triangles-1 sh -c "${BIN} domain.mesh")
add_test(convert_to_triangles-2 sh -c "${BIN} bracket.mesh")
add_test(convert_to_triangles-3 sh -c "${BIN} square.mesh")
add_test(convert_to_triangles-4 sh -c "${BIN} square_tri.mesh")

