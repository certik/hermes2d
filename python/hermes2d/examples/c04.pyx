from hermes2d._hermes2d cimport scalar, H1Space, BC_ESSENTIAL

cdef int bc_type_04(int marker):
    return BC_ESSENTIAL

cdef scalar bc_values_04(int marker, double x, double y):
    return x*x + y*y

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_type_04)
    space.thisptr.set_bc_values(&bc_values_04)
