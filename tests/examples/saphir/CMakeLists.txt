if(NOT REAL)
    return()
endif(NOT REAL)
project(example-saphir)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-saphire ${BIN})
set_tests_properties(example-saphire PROPERTIES LABELS slow)
