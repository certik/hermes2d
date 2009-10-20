from hermes2d._hermes2d cimport scalar, FuncReal, GeomReal, WeakForm, \
        int_grad_u_grad_v, int_v, malloc, ExtDataReal, c_Ord, create_Ord, \
        FuncOrd, GeomOrd, ExtDataOrd

cdef scalar bilinear_form(int n, double *wt, FuncReal *u, FuncReal *v, GeomReal
        *e, ExtDataReal *ext):
    return int_grad_u_grad_v(n, wt, u, v)

cdef scalar linear_form_p2(int n, double *wt, FuncReal *u, GeomReal
        *e, ExtDataReal *ext):
    return 2*int_v(n, wt, u)

cdef scalar linear_form_m1(int n, double *wt, FuncReal *u, GeomReal
        *e, ExtDataReal *ext):
    return -int_v(n, wt, u)

cdef scalar linear_form_m4(int n, double *wt, FuncReal *u, GeomReal
        *e, ExtDataReal *ext):
    return -4*int_v(n, wt, u)

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

def set_forms(WeakForm wf, int k=2):
    wf.thisptr.add_biform(0, 0, &bilinear_form, &_order_bf)
    if k == 2:
        wf.thisptr.add_liform(0, &linear_form_p2, &_order_lf)
    elif k == -1:
        wf.thisptr.add_liform(0, &linear_form_m1, &_order_lf)
    elif k == -4:
        wf.thisptr.add_liform(0, &linear_form_m4, &_order_lf)
    else:
        raise NotImplementedError()
