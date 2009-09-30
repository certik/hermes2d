set(SRC
       hash.cpp mesh1.cpp mesh2.cpp regul.cpp refmap.cpp curved.cpp
       transform.cpp traverse.cpp
       shapeset.cpp precalc.cpp solution.cpp filter.cpp
       space.cpp space_h1.cpp space_hcurl.cpp space_l2.cpp
       space_hdiv.cpp
       linear1.cpp linear2.cpp linear3.cpp graph.cpp
       view1.cpp view2.cpp view3.cpp view4.cpp view5.cpp
       quad_std.cpp
       shapeset_h1_ortho.cpp shapeset_h1_beuchler.cpp shapeset_h1_quad.cpp
       shapeset_hc_legendre.cpp shapeset_hc_gradleg.cpp
       shapeset_hd_legendre.cpp
       shapeset_l2_legendre.cpp
       qsort.cpp norm.cpp
       adapt_ortho_h1.cpp adapt_ortho_hcurl.cpp adapt_ortho_l2.cpp
       common.cpp matrix.cpp hermes2d.cpp weakform.cpp linsystem.cpp
       nonlinsystem.cpp refsystem.cpp
       mesh_parser.cpp mesh_lexer.cpp

       compat/fmemopen.cpp
       )

find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})

if(WITH_GLUT)
    find_package(GLUT REQUIRED)
endif(WITH_GLUT)

macro(LIBRARY_SETTINGS BIN)
    target_link_libraries(${BIN} ${JUDY_LIBRARY} pthread)
    if(WITH_GLUT)
        target_link_libraries(${BIN} glut)
    endif(WITH_GLUT)
endmacro(LIBRARY_SETTINGS)


INCLUDE (CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES m)
CHECK_FUNCTION_EXISTS(fmemopen HAVE_FMEMOPEN)
CHECK_FUNCTION_EXISTS(log2 HAVE_LOG2)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.h)


# build real version:
if(REAL)
    if(DEBUG)
        add_library(hermes2d-real-debug SHARED ${SRC})
        set_target_properties(hermes2d-real-debug PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
        LIBRARY_SETTINGS(hermes2d-real-debug)
        install(TARGETS hermes2d-real-debug LIBRARY DESTINATION lib)
    endif(DEBUG)
    if(RELEASE)
        add_library(hermes2d-real SHARED ${SRC})
        set_target_properties(hermes2d-real PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
        LIBRARY_SETTINGS(hermes2d-real)
        install(TARGETS hermes2d-real
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib
        )
    endif(RELEASE)
endif(REAL)

# build complex version:
if(COMPLEX)
    if(DEBUG)
        add_library(hermes2d-cplx-debug SHARED ${SRC})
        set_target_properties(hermes2d-cplx-debug PROPERTIES
            COMPILE_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -DCOMPLEX")
        LIBRARY_SETTINGS(hermes2d-cplx-debug)
        install(TARGETS hermes2d-cplx-debug LIBRARY DESTINATION lib)
    endif(DEBUG)
    if(RELEASE)
        add_library(hermes2d-cplx SHARED ${SRC})
        set_target_properties(hermes2d-cplx PROPERTIES
            COMPILE_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -DCOMPLEX")
        LIBRARY_SETTINGS(hermes2d-cplx)
        install(TARGETS hermes2d-cplx LIBRARY DESTINATION lib)
    endif(RELEASE)
endif(COMPLEX)

add_custom_target(debug)
add_dependencies(debug
    hermes2d-real-debug
    hermes2d-cplx-debug)

add_custom_target(release)
add_dependencies(release
    hermes2d-real
    hermes2d-cplx)

file(GLOB header_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${header_files} DESTINATION include/hermes2d)

