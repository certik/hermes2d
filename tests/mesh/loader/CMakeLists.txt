project(mesh-loader)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

add_test(h2d-1 sh -c "${BIN} h2d domain.mesh | diff - domain.dump")

if (WITH_EXODUSII)
  add_test(exodus-ii-1 sh -c "${BIN} exII square.e | diff - square.dump")
  add_test(exodus-ii-2 sh -c "${BIN} exII square-tri.e | diff - square-tri.dump")
  add_test(exodus-ii-3 sh -c "${BIN} exII iron-water.e | diff - iron-water.dump")
endif(WITH_EXODUSII)
