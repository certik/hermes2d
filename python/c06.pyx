from hermes2d cimport scalar, RealFunction, RefMap, DiscreteProblem, \
        int_grad_u_grad_v, int_v, H1Space, EdgePos, surf_int_G_v, surf_int_v, \
        surf_int_u_v, BC_ESSENTIAL, BC_NATURAL

cdef scalar bilinear_form(RealFunction *fu, RealFunction *fv,
        RefMap *ru, RefMap *rv):
    return int_grad_u_grad_v(fu, fv, ru, rv)

cdef scalar bilinear_form_surf_06(RealFunction *fu, RealFunction *fv,
        RefMap *ru, RefMap *rv, EdgePos* ep):
    if ep.marker != 1:
        return 0.
    return surf_int_u_v(fu, fv, ru, rv, ep)

cdef scalar linear_form_surf_06(RealFunction *fv, RefMap *rv, EdgePos* ep):
    if ep.marker != 1:
        return 0.
    return 20. * surf_int_v(fv, rv, ep)

cdef int bc_type_06(int marker):
    if marker == 3:
        return BC_ESSENTIAL
    else:
        return BC_NATURAL

cdef scalar bc_values_06(int marker, double x, double y):
    if marker == 3:
        return 100.
    return 0.

def set_forms(DiscreteProblem dp):
    dp.thisptr.set_bilinear_form(0, 0, &bilinear_form, NULL,
            &bilinear_form_surf_06)
    dp.thisptr.set_linear_form(0, NULL, &linear_form_surf_06);

def set_bc(H1Space space):
    space.thisptr.set_bc_types(&bc_type_06)
    space.thisptr.set_bc_values(&bc_values_06)
