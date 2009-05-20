from hermes2d.hermes2d cimport scalar, RealFunction, RefMap, WeakForm, \
        int_grad_u_grad_v, int_v, malloc

cdef struct Data:
    int k

cdef scalar bilinear_form(RealFunction *fu, RealFunction *fv,
        RefMap *ru, RefMap *rv):
    return int_grad_u_grad_v(fu, fv, ru, rv)

cdef scalar linear_form(RealFunction *fv, RefMap *rv, void *data):
    cdef Data *d = <Data *>data
    return d.k * int_v(fv, rv)

def set_forms(WeakForm wf, int k=2):
    wf.thisptr.add_biform(0, 0, &bilinear_form)
    wf.thisptr.add_liform(0, &linear_form)
    cdef Data *data
    data = <Data *>malloc(sizeof(Data))
    data.k = k
    wf.thisptr.add_liform_data(0, data)
