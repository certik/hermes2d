project(zoom-to-fit)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(zoom-to-fit_const-on-square "${BIN}" 0 square.mesh)
add_test(zoom-to-fit_const-on-distorted "${BIN}" 0 distorted.mesh)
add_test(zoom-to-fit_plane-on-square "${BIN}" 1 square.mesh)
add_test(zoom-to-fit_plane-on-distorted "${BIN}" 1 distorted.mesh)
add_test(zoom-to-fit_cuboid-on-square "${BIN}" 2 square.mesh)
add_test(zoom-to-fit_cuboid-on-curved "${BIN}" 2 curved.mesh)
add_test(zoom-to-fit_bigrange-on-square "${BIN}" 3 square.mesh)
add_test(zoom-to-fit_bigrange-on-curved "${BIN}" 3 curved.mesh)
add_test(zoom-to-fit_paraboloid-on-distorted "${BIN}" 4 distorted.mesh)

