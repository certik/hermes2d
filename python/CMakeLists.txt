include(UseCython)

find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
include_directories(${JUDY_INCLUDE_DIR})
include_directories(../src)

add_subdirectory(hermes2d)

add_custom_target(test
    COMMAND ./test tests
    )
