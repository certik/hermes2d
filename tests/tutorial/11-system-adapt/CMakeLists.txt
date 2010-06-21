project(tutorial-11-system-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-11-system-adapt-projection ${BIN})
set_tests_properties(tutorial-11-system-adapt-projection PROPERTIES LABELS slow)
add_test(tutorial-11-system-adapt-coarse ${BIN} -coarse_mesh)
set_tests_properties(tutorial-11-system-adapt-coarse PROPERTIES LABELS slow)
