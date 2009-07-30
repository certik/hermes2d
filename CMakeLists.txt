project(hermes2d)

cmake_minimum_required(VERSION 2.6)

#---------------------------------------------------------
# begin of configuration section
# The values below are sane defaults to produce a build on all platforms,
# if you want to override them, set them in CMake.vars


# build with debugging symbols
set(DEBUG NO)
set(RELEASE YES)

# real/complex version of the library (you can build either of both)
set(REAL     YES)
set(COMPLEX  YES)

# optional functionality
set(WITH_PYTHON   NO)
set(WITH_EXAMPLES YES)
set(WITH_TUTORIAL YES)
set(WITH_TESTS    NO)
set(WITH_GLUT     YES)
set(WITH_UTIL     YES)

# allow to override the default values in CMake.vars
if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
	include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)

# end of configuration section
#---------------------------------------------------------


if(REAL)
    list(APPEND VERS "real")
endif(REAL)
if(COMPLEX)
    list(APPEND VERS "complex")
endif(COMPLEX)

if(DEBUG)
    if(REAL)
        set(HERMES_REAL_BIN hermes2d-real-debug)
    endif(REAL)
    if(COMPLEX)
        set(HERMES_CPLX_BIN hermes2d-cplx-debug)
    endif(COMPLEX)
endif(DEBUG)

if(RELEASE)
    if(REAL)
        set(HERMES_REAL_BIN hermes2d-real)
    endif(REAL)
    if(COMPLEX)
        set(HERMES_CPLX_BIN hermes2d-cplx)
    endif(COMPLEX)
endif(RELEASE)

if(NOT WITH_GLUT)
    add_definitions(-DNOGLUT)
endif(NOT WITH_GLUT)


set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

add_subdirectory(src)

if(WITH_EXAMPLES)
	add_subdirectory(examples)
endif(WITH_EXAMPLES)

if(WITH_TUTORIAL)
    add_subdirectory(tutorial EXCLUDE_FROM_ALL)
    add_custom_target(tutorial)
    add_custom_command(
        TARGET	tutorial
        COMMAND make -C ${CMAKE_BINARY_DIR}/tutorial
    )
endif(WITH_TUTORIAL)

if(WITH_PYTHON)
	add_subdirectory(python)
endif(WITH_PYTHON)

if(WITH_UTIL)
	add_subdirectory(util)
endif(WITH_UTIL)

# Documentation
# This doesn't work yet:
#add_subdirectory(doc)
#add_custom_target(doc)
#add_custom_command(
#	SOURCE	${DOXYGEN_CONFIG_FILE}
#	COMMAND	${DOXYGEN_BINARY}
#	ARGS	${DOXYGEN_CONFIG_FILE}
#	TARGET	doc
#	OUTPUTS	${PROJECT_BINARY_DIR}/doc/html
#)

#add_custom_command(
#	SOURCE	doc
#	TARGET	doc
#	DEPENDS	${PROJECT_BINARY_DIR}/doc/html
#)


message("\nConfiguration results (these options can be set on the commandline)")
message("---------------------")
message("Versions to build: ${VERS}")
message("Build with Python: ${WITH_PYTHON}")
message("Build with examples: ${WITH_EXAMPLES}")
message("Build with tutorial: ${WITH_TUTORIAL}")
message("Build with GLUT: ${WITH_GLUT}")
message("Build with util: ${WITH_UTIL}")
message("\n")
