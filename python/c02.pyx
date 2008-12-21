from hermes2d cimport scalar, RealFunction, RefMap, DiscreteProblem, \
        int_grad_u_grad_v, int_v

cdef scalar bilinear_form(RealFunction *fu, RealFunction *fv,
        RefMap *ru, RefMap *rv):
    return int_grad_u_grad_v(fu, fv, ru, rv)

cdef scalar linear_form(RealFunction *fv, RefMap *rv):
    return 2*int_v(fv, rv)

def set_forms(DiscreteProblem dp):
    dp.thisptr.set_bilinear_form(0, 0, &bilinear_form)
    dp.thisptr.set_linear_form(0, &linear_form);
