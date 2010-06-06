project(tutorial-10-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-10-adapt-projection ${BIN})
set_tests_properties(tutorial-10-adapt-projection PROPERTIES LABELS slow)
add_test(tutorial-10-adapt-coarse ${BIN} -coarse_mesh)
set_tests_properties(tutorial-10-adapt-coarse PROPERTIES LABELS slow)
