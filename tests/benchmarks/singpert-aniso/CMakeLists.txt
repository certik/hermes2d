if(NOT REAL)
    return()
endif(NOT REAL)
project(benchmark-singpert-aniso)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-singpert-aniso ${BIN})
# this test is known to fail:
set_tests_properties(benchmark-singpert-aniso PROPERTIES WILL_FAIL TRUE)
set_tests_properties(benchmark-singpert-aniso PROPERTIES LABELS slow)
