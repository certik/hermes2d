from hermes2d._hermes2d cimport scalar, H1Space, BC_ESSENTIAL

# Boundary condition type (essential = Dirichlet)
cdef int bc_type_04(int marker):
    return BC_ESSENTIAL

# Function values for Dirichlet boundary markers
cdef scalar bc_values_04(int marker, double x, double y):
    return x*x + y*y

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_type_04)
    space.thisptr.set_bc_values(&bc_values_04)

"""
cdef scalar bilinear_form(int n, double *wt, FuncReal *u, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return result

# Volumetric linear form (right-hand side)
cdef scalar linear_form(int n, double *wt, FuncReal *v, GeomReal *e, ExtDataReal *ext):
    return int_F_v(n, wt, rhs, v, e)


# Integration order for the volumetric linear form
cdef c_Ord linear_form_ord(int n, double *wt, FuncOrd *v, GeomOrd *e, ExtDataOrd *ext):
    #return v.val[0] * e.x[0] * e.x[0]
    return create_Ord(20)

def set_forms(WeakForm wf):
    wf.thisptr.add_biform(0, 0, &bilinear_form, &bilinear_form_ord, H2D_SYM);
    wf.thisptr.add_liform(0, &linear_form, &linear_form_ord);
    wf.thisptr.add_liform_surf(0, &linear_form_surf, &linear_form_surf_ord, 2);
"""
