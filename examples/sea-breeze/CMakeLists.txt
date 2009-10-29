if(NOT REAL)
    return()
endif(NOT REAL)

project(sea-breeze)

find_package(PythonLibs REQUIRED)

include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
include_directories(${hermes2d_SOURCE_DIR}/python/hermes2d)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)


target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES})
