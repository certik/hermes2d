add_subdirectory(examples)
cython_add_module(_hermes2d dummy_solver.cpp stdcython.c)
target_link_libraries(_hermes2d ${HERMES_REAL_BIN})
if(WITH_PYTHON_UMFPACK)
    find_package(UMFPACK REQUIRED)
    include_directories(${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})
    target_link_libraries(_hermes2d umfpack amd blas)
    set(USE_UMFPACK 1)
else(WITH_PYTHON_UMFPACK)
    set(USE_UMFPACK 0)
endif(WITH_PYTHON_UMFPACK)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

cython_add_module(forms)
target_link_libraries(forms ${HERMES_REAL_BIN})
install(TARGETS forms LIBRARY DESTINATION lib/python/hermes2d)
set_target_properties(forms PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/python/hermes2d"
    )

install(TARGETS _hermes2d LIBRARY DESTINATION lib/python/hermes2d)
set_target_properties(_hermes2d PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/python/hermes2d"
    )

install(FILES __init__.py plot.py runtests.py mesh.py
    DESTINATION lib/python/hermes2d
    )

# tests
install(FILES tests/test_basics.py tests/test_examples.py tests/test_plots.py
    tests/test_mesh.py
    DESTINATION lib/python/hermes2d/tests
    )

install(FILES _hermes2d.pxd dummy_solver.h config.h stdcython.h
    DESTINATION lib/python/hermes2d/include
    )
