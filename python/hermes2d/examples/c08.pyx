from hermes2d._hermes2d cimport scalar, RealFunction, RefMap, WeakForm, \
        int_grad_u_grad_v, int_v, H1Space, Solution, int_u_dvdx, \
        int_u_dvdy, int_w_nabla_u_v, int_u_v, BF_ANTISYM, BC_ESSENTIAL, \
        BC_NONE, SYM, UNSYM, ANY, ANTISYM
from hermes2d._hermes2d cimport scalar, WeakForm, H1Space, EdgePos, \
        FuncReal, GeomReal, ExtDataReal, BC_ESSENTIAL, \
        BC_NATURAL, int_v, c_Ord, create_Ord, FuncOrd, GeomOrd, ExtDataOrd
from hermes2d._hermes2d cimport int_dudy_dvdy, int_dudy_dvdx, int_dudx_dvdx, \
        int_dudx_dvdy, SYM

cdef Solution xprev, yprev
cdef double Re = 700
cdef double tau = 0.05

cdef int xvel_bc_type(int marker):
    if marker != 2:
        return BC_ESSENTIAL
    else:
        return BC_NONE

cdef scalar xvel_bc_value(int marker, double x, double y):
    if marker != 5:
        return 1
    else:
        return 0

cdef int yvel_bc_type(int marker):
    if marker != 2:
        return BC_ESSENTIAL
    else:
        return BC_NONE

cdef int press_bc_type(int marker):
    return BC_NONE

cdef scalar bilinear_form_sym_0_0_1_1(int n, double *wt, FuncReal *u, FuncReal
        *v, GeomReal *e, ExtDataReal *ext):
    return int_grad_u_grad_v(n, wt, u, v)/Re + int_u_v(n, wt, u, v)/tau

cdef scalar bilinear_form_unsym_0_0_1_1(int n, double *wt, FuncReal *u, FuncReal
        *v, GeomReal *e, ExtDataReal *ext):
    return int_w_nabla_u_v(n, wt, ext.fn[0], ext.fn[1], u, v)

cdef scalar bilinear_form_unsym_0_2(int n, double *wt, FuncReal *u, FuncReal
        *v, GeomReal *e, ExtDataReal *ext):
    return -int_u_dvdx(n, wt, u, v)

cdef scalar bilinear_form_unsym_1_2(int n, double *wt, FuncReal *u, FuncReal
        *v, GeomReal *e, ExtDataReal *ext):
    return -int_u_dvdy(n, wt, u, v)

cdef scalar linear_form(int n, double *wt, FuncReal *v, GeomReal *e, ExtDataReal
        *ext):
    return int_u_v(n, wt, ext.fn[0], v)/tau

cdef c_Ord _order_bf(int n, double *wt, FuncOrd *u, FuncOrd *v, GeomOrd
        *e, ExtDataOrd *ext):
    # XXX: with 9 it doesn't shout about the integration order, but gives wrong
    # results...
    return create_Ord(20)

cdef c_Ord _order_lf(int n, double *wt, FuncOrd *u, GeomOrd
        *e, ExtDataOrd *ext):
    # XXX: with 9 it doesn't shout about the integration order, but gives wrong
    # results...
    return create_Ord(20)

def set_forms(WeakForm wf, Solution xprev2, Solution yprev2):
    global xprev
    xprev = xprev2
    global yprev
    yprev = yprev2
    wf.thisptr.add_biform(0, 0, &bilinear_form_sym_0_0_1_1, &_order_bf, SYM)
    wf.thisptr.add_biform(0, 0, &bilinear_form_unsym_0_0_1_1, &_order_bf,
            UNSYM, ANY, 2, xprev.thisptr, yprev.thisptr)
    wf.thisptr.add_biform(1, 1, &bilinear_form_sym_0_0_1_1, &_order_bf, SYM)
    wf.thisptr.add_biform(1, 1, &bilinear_form_unsym_0_0_1_1, &_order_bf,
            UNSYM, ANY, 2, xprev.thisptr, yprev.thisptr)
    wf.thisptr.add_biform(0, 2, &bilinear_form_unsym_0_2, &_order_bf, ANTISYM)
    wf.thisptr.add_biform(1, 2, &bilinear_form_unsym_1_2, &_order_bf, ANTISYM)
    wf.thisptr.add_liform(0, &linear_form, &_order_lf, ANY, 1, xprev.thisptr)
    wf.thisptr.add_liform(1, &linear_form, &_order_lf, ANY, 1, yprev.thisptr)

def set_bc(H1Space xvel, H1Space yvel, H1Space press):
    xvel.thisptr.set_bc_types(&xvel_bc_type)
    xvel.thisptr.set_bc_values(&xvel_bc_value)
    yvel.thisptr.set_bc_types(&yvel_bc_type)
    press.thisptr.set_bc_types(&press_bc_type)
