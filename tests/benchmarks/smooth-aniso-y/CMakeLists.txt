if(NOT REAL)
    return()
endif(NOT REAL)
project(benchmark-smooth-aniso-y)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-smooth-aniso-y ${BIN})
# this test is known to fail:
set_tests_properties(benchmark-smooth-aniso-y PROPERTIES WILL_FAIL TRUE)
