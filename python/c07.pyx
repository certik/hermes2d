from hermes2d cimport scalar, RealFunction, RefMap, DiscreteProblem, \
        int_grad_u_grad_v, int_v, H1Space, EdgePos, surf_int_G_v, surf_int_v, \
        surf_int_u_v, BF_SYM, int_a_dudx_dvdy_b_dudy_dvdx, \
        int_a_dudx_dvdx_b_dudy_dvdy, BC_ESSENTIAL, BC_NATURAL

cdef double E  = 200e9
cdef double nu = 0.3
cdef double l = (E * nu) / ((1 + nu) * (1 - 2*nu))
cdef double mu = E / (2*(1 + nu))

cdef scalar bilinear_form_07_0_0(RealFunction* fu, RealFunction* fv,
    RefMap* ru, RefMap* rv):
    return int_a_dudx_dvdx_b_dudy_dvdy(l+2*mu, fu, mu, fv, ru, rv)

cdef scalar bilinear_form_07_0_1(RealFunction* fu, RealFunction* fv,
    RefMap* ru, RefMap* rv):
    return int_a_dudx_dvdy_b_dudy_dvdx(l, fv, mu, fu, rv, ru)

cdef scalar bilinear_form_07_1_1(RealFunction* fu, RealFunction* fv,
    RefMap* ru, RefMap* rv):
    return int_a_dudx_dvdx_b_dudy_dvdy(mu, fu, l+2*mu, fv, ru, rv)

cdef scalar linear_form_surf_07(RealFunction* fv, RefMap* rv, EdgePos *ep):
    return surf_int_G_v(fv, rv, ep)

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

def set_forms(DiscreteProblem dp):
    dp.thisptr.set_bilinear_form(0, 0, NULL, &bilinear_form_07_0_0)
    dp.thisptr.set_bilinear_form(1, 1, NULL, &bilinear_form_07_1_1)
    dp.thisptr.set_bilinear_form(0, 1, &bilinear_form_07_0_1)
    dp.thisptr.set_bilinear_form(1, 0, BF_SYM)
    dp.thisptr.set_linear_form(1, NULL, &linear_form_surf_07);

def set_bc(H1Space xdisp, H1Space ydisp):
    xdisp.thisptr.set_bc_types(&bc_type_07x)
    ydisp.thisptr.set_bc_types(&bc_type_07y)
    ydisp.thisptr.set_bc_values_edge(&bc_values_07y)
