project(tutorial-11-adapt-system)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-11-adapt-system-projection ${BIN})
set_tests_properties(tutorial-11-adapt-system-projection PROPERTIES LABELS slow)
add_test(tutorial-11-adapt-system-coarse ${BIN} -coarse_mesh)
set_tests_properties(tutorial-11-adapt-system-coarse PROPERTIES LABELS slow)
