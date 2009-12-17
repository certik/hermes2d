project(convert_to_quads)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(convert_to_quads-1 sh -c "${BIN} domain.mesh")
add_test(convert_to_quads-2 sh -c "${BIN} bracket.mesh")
add_test(convert_to_quads-3 sh -c "${BIN} square.mesh")
add_test(convert_to_quads-4 sh -c "${BIN} square_tri.mesh")

