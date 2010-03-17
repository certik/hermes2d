from hermes2d._hermes2d cimport scalar, WeakForm, H1Space, EdgePos, \
        FuncReal, GeomReal, ExtDataReal, BC_ESSENTIAL, \
        BC_NATURAL, c_Ord, create_Ord, FuncOrd, GeomOrd, \
        ExtDataOrd, ExtDataReal, FuncReal, GeomReal, SYM, int_F_v

def set_bc_x(H1Space space):
    pass

def set_bc_y(H1Space space):
    pass

def set_bc(H1Space space):
    pass

def set_forms(WeakForm wf):
    pass

"""
# Boundary condition types for x-velocity
cdef int xvel_bc_type(int marker):
    if (marker == marker_right):
        return BC_NONE
    else:
        return BC_ESSENTIAL
        
# Boundary condition values for x-velocity
cdef scalar xvel_bc_value(int marker, double x, double y):
    if (marker == marker_left):
        # Time-dependent inlet velocity (parabolic profile)
        val_y = VEL_INLET * y*(H-y) / (H/2.)/(H/2.)
        if (TIME <= STARTUP_TIME):
            return val_y * TIME/STARTUP_TIME
        else:
            return val_y
    else:
        return 0

# Boundary condition types for y-velocity
cdef int yvel_bc_type(int marker):
    if (marker == marker_right):
        return BC_NONE
    else:
        return BC_ESSENTIAL

cdef int p_bc_type(int marker):
    return BC_NONE

def set_bc_x(H1Space space):
    space.thisptr.set_bc_types(&xvel_bc_type)
    space.thisptr.set_bc_values(&xvel_bc_type)

def set_bc_y(H1Space space):
    space.thisptr.set_bc_types(&yvel_bc_type)

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&p_bc_type)
"""
