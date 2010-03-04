IF(COMPLEX)

project(newton-timedep-gp)

# use the complex version of the library:
set(HERMES ${HERMES_CPLX_BIN})

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)

ENDIF(COMPLEX)
