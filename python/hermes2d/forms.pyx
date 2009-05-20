from hermes2d.hermes2d cimport scalar, RealFunction, RefMap, WeakForm, \
        int_grad_u_grad_v, int_v, malloc

cdef scalar bilinear_form(RealFunction *fu, RealFunction *fv,
        RefMap *ru, RefMap *rv):
    return int_grad_u_grad_v(fu, fv, ru, rv)

cdef scalar linear_form_p2(RealFunction *fv, RefMap *rv):
    return 2*int_v(fv, rv)

cdef scalar linear_form_m1(RealFunction *fv, RefMap *rv):
    return -int_v(fv, rv)

cdef scalar linear_form_m4(RealFunction *fv, RefMap *rv):
    return -4*int_v(fv, rv)

def set_forms(WeakForm wf, int k=2):
    wf.thisptr.add_biform(0, 0, &bilinear_form)
    if k == 2:
        wf.thisptr.add_liform(0, &linear_form_p2)
    elif k == -1:
        wf.thisptr.add_liform(0, &linear_form_m1)
    elif k == -4:
        wf.thisptr.add_liform(0, &linear_form_m4)
    else:
        raise NotImplementedError()
