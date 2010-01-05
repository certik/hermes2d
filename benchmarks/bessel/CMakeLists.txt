IF(COMPLEX)

project(bessel)

# use the complex version of the library:
set(HERMES ${HERMES_CPLX_BIN})

add_executable(${PROJECT_NAME} main.cpp bessel.cpp)
include (../CMake.common)

ENDIF(COMPLEX)
