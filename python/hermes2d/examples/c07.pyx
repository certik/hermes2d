from hermes2d._hermes2d cimport scalar, WeakForm, H1Space, EdgePos, \
        FuncReal, GeomReal, ExtDataReal, BC_ESSENTIAL, \
        BC_NATURAL, int_v, c_Ord, create_Ord, FuncOrd, GeomOrd, ExtDataOrd
from hermes2d._hermes2d cimport int_dudy_dvdy, int_dudy_dvdx, int_dudx_dvdx, \
        int_dudx_dvdy, SYM, int_v_ord

cdef double E  = 200e9
cdef double nu = 0.3
cdef double l = (E * nu) / ((1 + nu) * (1 - 2*nu))
cdef double mu = E / (2*(1 + nu))
cdef double f = 1e4

cdef scalar bilinear_form_0_0(int n, double *wt, FuncReal *u, FuncReal *v,
        GeomReal *e, ExtDataReal *ext):
    return (l +2*mu) * int_dudx_dvdx(n, wt, u, v) + \
            mu * int_dudy_dvdy(n, wt, u, v)

cdef scalar bilinear_form_0_1(int n, double *wt, FuncReal *u, FuncReal *v,
        GeomReal *e, ExtDataReal *ext):
    return l * int_dudy_dvdx(n, wt, u, v) + mu * int_dudx_dvdy(n, wt, u, v)

cdef scalar bilinear_form_1_1(int n, double *wt, FuncReal *u, FuncReal *v,
        GeomReal *e, ExtDataReal *ext):
    return mu * int_dudx_dvdx(n, wt, u, v) + \
            (l+2*mu) * int_dudy_dvdy(n, wt, u, v)

cdef scalar linear_form_surf_07(int n, double *wt, FuncReal *u, GeomReal
        *e, ExtDataReal *ext):
    return f * int_v(n, wt, u)

cdef int bc_type_07x(int marker):
    return BC_NATURAL

cdef int bc_type_07y(int marker):
    if marker == 1:
        return BC_ESSENTIAL
    return BC_NATURAL

cdef scalar bc_values_07y(EdgePos *ep):
    if ep.marker == 3:
        return 1e4
    return 0.

cdef c_Ord _order_bf(int n, double *wt, FuncOrd *u, FuncOrd *v, GeomOrd
        *e, ExtDataOrd *ext):
    # XXX: with 9 it doesn't shout about the integration order, but gives wrong
    # results...
    return create_Ord(20)

cdef c_Ord _order_lf(int n, double *wt, FuncOrd *u, GeomOrd
        *e, ExtDataOrd *ext):
    return int_v_ord(n, wt, u).mul_double(f)

def set_forms(WeakForm dp):
    dp.thisptr.add_biform(0, 0, &bilinear_form_0_0, &_order_bf, SYM)
    dp.thisptr.add_biform(0, 1, &bilinear_form_0_1, &_order_bf, SYM)
    dp.thisptr.add_biform(1, 1, &bilinear_form_1_1, &_order_bf, SYM)
    dp.thisptr.add_liform_surf(1, &linear_form_surf_07, &_order_lf);

def set_bc(H1Space xdisp, H1Space ydisp):
    xdisp.thisptr.set_bc_types(&bc_type_07x)
    ydisp.thisptr.set_bc_types(&bc_type_07y)
    ydisp.thisptr.set_bc_values_edge(&bc_values_07y)
