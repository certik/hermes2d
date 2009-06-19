set(SRC
       hash.cpp mesh1.cpp mesh2.cpp regul.cpp refmap.cpp curved.cpp
       transform.cpp traverse.cpp
       shapeset.cpp precalc.cpp solution.cpp filter.cpp
       space.cpp space_h1.cpp space_hcurl.cpp space_l2.cpp
       space_hdiv.cpp
       linear1.cpp linear2.cpp linear3.cpp graph.cpp
       view1.cpp view2.cpp view3.cpp view4.cpp
       quad_std.cpp
       shapeset_h1_ortho.cpp shapeset_h1_beuchler.cpp shapeset_h1_quad.cpp
       shapeset_hc_legendre.cpp shapeset_hc_gradleg.cpp
       shapeset_hd_legendre.cpp
       qsort.cpp norm.cpp
       adapt_ortho_h1.cpp adapt_ortho_hcurl.cpp adapt_ortho_l2.cpp
       common.cpp matrix.cpp hermes2d.cpp weakform.cpp linsystem.cpp
       nonlinsystem.cpp refsystem.cpp
       mesh_parser.cpp mesh_lexer.cpp
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


# build real version:
if(REAL)
    add_library(${HERMES_REAL_BIN} SHARED ${SRC})
    LIBRARY_SETTINGS(${HERMES_REAL_BIN})
endif(REAL)

# build complex version:
if(COMPLEX)
    add_library(${HERMES_CPLX_BIN} SHARED ${SRC})
    set_target_properties(${HERMES_CPLX_BIN} PROPERTIES COMPILE_FLAGS -DCOMPLEX)
    LIBRARY_SETTINGS(${HERMES_CPLX_BIN})
endif(COMPLEX)

install(TARGETS ${HERMES_REAL_BIN} LIBRARY DESTINATION lib)
file(GLOB header_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${header_files} DESTINATION include/hermes2d)
