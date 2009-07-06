from hermes2d._hermes2d cimport scalar, RealFunction, RefMap, WeakForm, \
        H1Space, EdgePos, surf_int_G_v, \
        BC_ESSENTIAL, BC_NATURAL, c_atan, c_pi, c_sqrt, c_sqr, int_F_v, \
        int_grad_u_grad_v

cdef scalar linear_form_surf_05(RealFunction *fv, RefMap *rv, EdgePos* ep):
    return surf_int_G_v(fv, rv, ep)

cdef double fn(double x, double y):
    return c_atan(60 * (c_sqrt(c_sqr(x-1.25) + c_sqr(y+0.25)) - c_pi/3))

cdef double rhs(double x, double y):
    cdef double t1 = c_sqrt(16*(x*x + y*y) - 40*x + 8*y + 26)
    cdef double t2 = 3600*(x*x + y*y) - 9000*x + 1800*y
    return -(240 * (t2 + 5849 - 400*c_pi*c_pi)) / (t1 * c_sqr(5851 + t2 - \
              600*t1*c_pi + 400*c_pi*c_pi))

cdef scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru,
        RefMap* rv):
    return int_grad_u_grad_v(fu, fv, ru, rv)

cdef scalar linear_form(RealFunction *fv, RefMap *rv):
    return -int_F_v(rhs, fv, rv)

cdef scalar bc_values(int marker, double x, double y):
    return fn(x, y)

def set_forms(WeakForm dp):
    dp.thisptr.add_biform(0, 0, &bilinear_form)
    dp.thisptr.add_liform(0, &linear_form)

def set_bc(H1Space space):
    space.thisptr.set_bc_values(&bc_values)
