from hermes2d._hermes2d cimport scalar, H1Space, BC_ESSENTIAL

# boundary condition types (essential = Dirichlet)
cdef int bc_type(int marker):
    return BC_ESSENTIAL

# function values for Dirichlet boundary conditions
cdef scalar bc_values(int marker, double x, double y):
    return 0

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_type)
    space.thisptr.set_bc_values(&bc_values)
