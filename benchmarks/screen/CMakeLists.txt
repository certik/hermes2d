IF(COMPLEX)

project(screen)

# use the complex version of the library:
set(HERMES ${HERMES_CPLX_BIN})

add_executable(${PROJECT_NAME} main.cpp fresnl.c)
include (../CMake.common)

ENDIF(COMPLEX)
