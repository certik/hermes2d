project(hermes_common)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

include(UseCython)
# if the user specified PYTHON_INCLUDE_PATH, let's use that, otherwise let's
# try to guess using cmake/FindPythonLibs.cmake
if (NOT PYTHON_INCLUDE_PATH)
    find_package(PythonLibs REQUIRED)
endif (NOT PYTHON_INCLUDE_PATH)
find_package(NumPy REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})

enable_testing()

#PYTHONPATH=${hermes_common_SOURCE_DIR}
set(HERMES_COMMON hermes_common)





set(SRC
    _hermes_common_api_new.cpp
    matrix.cpp
    python_solvers.cpp
    python_api.cpp
    )

add_library(${HERMES_COMMON} SHARED ${SRC})
target_link_libraries(${HERMES_COMMON} ${PYTHON_LIBRARIES})
set_source_files_properties(matrix.cpp PROPERTIES
    OBJECT_DEPENDS ${hermes_common_SOURCE_DIR}/_hermes_common_api_new.h
    )
add_custom_command(
    OUTPUT _hermes_common_api.h
    COMMAND cython _hermes_common.pyx
    )
add_subdirectory(tests)

cython_add_module(_hermes_common stdcython.cpp)
target_link_libraries(_hermes_common ${HERMES_COMMON})
set_property(TARGET _hermes_common PROPERTY COMPILE_FLAGS "-g")
