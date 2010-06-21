project(tutorial-12-general-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-12-general-adapt-projection ${BIN})
add_test(tutorial-12-general-adapt-coarse ${BIN} -coarse_mesh)

