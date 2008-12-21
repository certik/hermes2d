from hermes2d cimport scalar, RealFunction, RefMap, DiscreteProblem, \
        int_grad_u_grad_v, int_v, H1Space, Solution, int_u_dvdx, \
        int_u_dvdy, int_w_nabla_u_v, int_u_v, BF_ANTISYM, BC_ESSENTIAL, \
        BC_NONE

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

cdef scalar bilinear_form_unsym_0_0_1_1(RealFunction* fu, RealFunction* fv,
    RefMap* ru, RefMap* rv):
    return int_w_nabla_u_v(<RealFunction *>(xprev.thisptr),
            <RealFunction *>(yprev.thisptr), fu, fv, ru, rv)

cdef scalar bilinear_form_unsym_0_2(RealFunction* fu, RealFunction* fv,
    RefMap* ru, RefMap* rv):
    return -int_u_dvdx(fu, fv, ru, rv)

cdef scalar bilinear_form_unsym_1_2(RealFunction* fu, RealFunction* fv,
    RefMap* ru, RefMap* rv):
    return -int_u_dvdy(fu, fv, ru, rv)

cdef scalar bilinear_form_sym_0_0_1_1(RealFunction* fu, RealFunction* fv,
    RefMap* ru, RefMap* rv):
  return int_grad_u_grad_v(fu, fv, ru, rv) / Re + int_u_v(fu, fv, ru, rv) / tau

cdef scalar linear_form_0(RealFunction* fv, RefMap* rv):
    return int_u_v(<RealFunction *>(xprev.thisptr), fv, rv, rv) / tau

cdef scalar linear_form_1(RealFunction* fv, RefMap* rv):
    return int_u_v(<RealFunction *>(yprev.thisptr), fv, rv, rv) / tau

def set_forms(DiscreteProblem dp, Solution xprev2, Solution yprev2):
    global xprev
    xprev = xprev2
    global yprev
    yprev = yprev2
    dp.thisptr.set_bilinear_form(0, 0, &bilinear_form_unsym_0_0_1_1,
            &bilinear_form_sym_0_0_1_1);
    dp.thisptr.set_bilinear_form(0, 2, &bilinear_form_unsym_0_2);
    dp.thisptr.set_bilinear_form(1, 1, &bilinear_form_unsym_0_0_1_1,
            &bilinear_form_sym_0_0_1_1);
    dp.thisptr.set_bilinear_form(1, 2, &bilinear_form_unsym_1_2);
    dp.thisptr.set_bilinear_form(2, 0, BF_ANTISYM);
    dp.thisptr.set_bilinear_form(2, 1, BF_ANTISYM);
    dp.thisptr.set_linear_form(0, &linear_form_0);
    dp.thisptr.set_linear_form(1, &linear_form_1);

def set_bc(H1Space xvel, H1Space yvel, H1Space press):
    xvel.thisptr.set_bc_types(&xvel_bc_type)
    xvel.thisptr.set_bc_values(&xvel_bc_value)
    yvel.thisptr.set_bc_types(&yvel_bc_type)
    press.thisptr.set_bc_types(&press_bc_type)
