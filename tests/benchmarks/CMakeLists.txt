find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})

# examples
add_subdirectory(nist-1) 
add_subdirectory(layer)
add_subdirectory(lshape)
add_subdirectory(smooth-iso)
add_subdirectory(smooth-aniso-x)
add_subdirectory(smooth-aniso-y)

#if(NOT WITH_TRILINOS)
  add_subdirectory(bessel)
  add_subdirectory(screen)
#endif(NOT WITH_TRILINOS)
