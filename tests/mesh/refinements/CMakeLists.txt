project(refinements)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(refinements-1 sh -c "${BIN} domain.mesh")
add_test(refinements-2 sh -c "${BIN} bracket.mesh")
add_test(refinements-3 sh -c "${BIN} square.mesh")
add_test(refinements-4 sh -c "${BIN} square_tri.mesh")
