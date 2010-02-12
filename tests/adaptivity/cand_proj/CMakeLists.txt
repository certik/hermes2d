project(cand_proj)

add_executable(${PROJECT_NAME} main.cpp functions.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(candidate_projection ${BIN})
add_test(candidate_projection ${BIN} -tri)
