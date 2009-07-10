if(NOT WITH_GLUT)
    message(FATAL_ERROR "util/linview must be compiled with GLUT")
endif(NOT WITH_GLUT)

if(NOT REAL)
    message(STATUS "skipping util/linview (real version is not being built)")
    return()
endif(NOT REAL)

macro(LINVIEW_SETTINGS VIEW)
    add_executable(${PROJECT_NAME} main.cpp)
    include_directories(${hermes2d_SOURCE_DIR}/src)
    target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS
        -DWHATVIEW=${VIEW})
endmacro(LINVIEW_SETTINGS)

project(linview)
LINVIEW_SETTINGS(ScalarView)

project(vecview)
LINVIEW_SETTINGS(VectorView)

project(ordview)
LINVIEW_SETTINGS(OrderView)
