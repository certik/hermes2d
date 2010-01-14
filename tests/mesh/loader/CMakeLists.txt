if (WITH_EXODUSII)
  project(mesh-loader)

  add_executable(${PROJECT_NAME} main.cpp)
  include (../../CMake.common)

  set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
  add_test(exodus-ii-1 sh -c "${BIN} square.e | diff - square.dump")
  add_test(exodus-ii-2 sh -c "${BIN} square-tri.e | diff - square-tri.dump")
  add_test(exodus-ii-3 sh -c "${BIN} iron-water.e | diff - iron-water.dump")
endif(WITH_EXODUSII)
