from hermes2d cimport scalar, RealFunction, RefMap, DiscreteProblem, \
        int_grad_u_grad_v, int_v, H1Space, EdgePos, surf_int_G_v, \
        BC_ESSENTIAL, BC_NATURAL

cdef scalar bilinear_form(RealFunction *fu, RealFunction *fv,
        RefMap *ru, RefMap *rv):
    return int_grad_u_grad_v(fu, fv, ru, rv)

cdef scalar linear_form_05(RealFunction *fv, RefMap *rv):
    return -int_v(fv, rv)

cdef scalar linear_form_surf_05(RealFunction *fv, RefMap *rv, EdgePos* ep):
    return surf_int_G_v(fv, rv, ep)

cdef int bc_type_05(int marker):
    if marker == 4:
        return BC_ESSENTIAL
    else:
        return BC_NATURAL

cdef scalar bc_values_05(int marker, double x, double y):
    if marker == 4:
        return 0.0
    elif marker == 2:
        return 1.0
    return 0.0

def set_forms(DiscreteProblem dp):
    dp.thisptr.set_bilinear_form(0, 0, &bilinear_form)
    dp.thisptr.set_linear_form(0, &linear_form_05, &linear_form_surf_05);

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_type_05)
    space.thisptr.set_bc_values(&bc_values_05)
