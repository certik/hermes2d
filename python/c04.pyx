from hermes2d cimport scalar, RealFunction, RefMap, DiscreteProblem, \
        int_grad_u_grad_v, int_v, H1Space, BC_ESSENTIAL

cdef scalar bilinear_form(RealFunction *fu, RealFunction *fv,
        RefMap *ru, RefMap *rv):
    return int_grad_u_grad_v(fu, fv, ru, rv)

cdef scalar linear_form(RealFunction *fv, RefMap *rv):
    return -4*int_v(fv, rv)

cdef int bc_type_04(int marker):
    return BC_ESSENTIAL

cdef scalar bc_values_04(int marker, double x, double y):
    return x*x + y*y

def set_forms(DiscreteProblem dp):
    dp.thisptr.set_bilinear_form(0, 0, &bilinear_form)
    dp.thisptr.set_linear_form(0, &linear_form);

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_type_04)
    space.thisptr.set_bc_values(&bc_values_04)
