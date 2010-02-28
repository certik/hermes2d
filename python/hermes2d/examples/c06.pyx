from hermes2d._hermes2d cimport scalar, FuncReal, GeomReal, WeakForm, \
        int_grad_u_grad_v, int_grad_u_grad_v_ord, int_v, int_v_ord, malloc, ExtDataReal, c_Ord, create_Ord, \
        FuncOrd, GeomOrd, ExtDataOrd, int_u_v, int_u_v_ord, BC_NATURAL, BC_ESSENTIAL, \
        H1Space

# Boundary condition types
cdef int bc_type_06(int marker):
    if marker == 3:
        return BC_ESSENTIAL
    else:
        return BC_NATURAL

# Function values for Dirichlet boundary markers
cdef scalar bc_values_06(int marker, double x, double y):
    if marker == 3:
        return 100.
    return 0.

cdef scalar bilinear_form(int n, double *wt, FuncReal *u, FuncReal *v, GeomReal
        *e, ExtDataReal *ext):
    return int_grad_u_grad_v(n, wt, u, v)

cdef scalar bilinear_form_surf_06(int n, double *wt, FuncReal *u, FuncReal *v,
        GeomReal *e, ExtDataReal *ext):
    if e.marker != 1:
        return 0.
    return int_u_v(n, wt, u, v)

cdef scalar linear_form_surf_06(int n, double *wt, FuncReal *v, GeomReal *e,
        ExtDataReal *ext):
    if e.marker != 1:
        return 0.
    return 20. * int_v(n, wt, v)

cdef c_Ord _order_bf(int n, double *wt, FuncOrd *u, FuncOrd *v, GeomOrd
        *e, ExtDataOrd *ext):
    return int_grad_u_grad_v_ord(n, wt, u, v)

cdef c_Ord _order_bf_surf_06(int n, double *wt, FuncOrd *u, FuncOrd *v,
        GeomOrd *e, ExtDataOrd *ext):
    return int_u_v_ord(n, wt, u, v)

cdef c_Ord _order_lf(int n, double *wt, FuncOrd *u, GeomOrd
        *e, ExtDataOrd *ext):
    return int_v_ord(n, wt, u)

def set_forms(WeakForm dp):
    dp.thisptr.add_biform(0, 0, &bilinear_form, &_order_bf)
    dp.thisptr.add_biform_surf(0, 0, &bilinear_form_surf_06, &_order_bf_surf_06)
    dp.thisptr.add_liform_surf(0, &linear_form_surf_06, &_order_lf)

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_type_06)
    space.thisptr.set_bc_values(&bc_values_06)
