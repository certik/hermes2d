from hermes2d._hermes2d cimport H1Space, c_H1Space, c_WeakForm, c_Solution, \
        WeakForm, Solution, c_Mesh, Mesh
from _hermes_common cimport c2numpy_double, numpy2c_double_inplace


cdef extern from "numerical_flux.h":
    double c_matrix_R "matrix_R"(int i, int j, double w0, double w1, double w3, double w4)
    double c_matrix_R_inv "matrix_R_inv"(int i, int j, double w0, double w1, double w3, double w4)
    double c_matrix_D_minus "matrix_D_minus"(int i, int j, double w0, double w1, double w3, double w4)
    void c_flux_riemann "riemann_solver"(double result[4], double w_l[4], double w_r[4])
    void c_flux_riemann_invert "riemann_solver_invert"(double result[4], double w_l[4], double w_r[4])
    void c_numerical_flux "numerical_flux"(double result[4], double w_l[4],
            double w_r[4], double nx, double ny)

def matrix_R(int i, int j, double w0, double w1, double w3, double w4):
    return c_matrix_R(i, j, w0, w1, w3, w4)

def matrix_R_inv(int i, int j, double w0, double w1, double w3, double w4):
    return c_matrix_R_inv(i, j, w0, w1, w3, w4)

def matrix_D_minus(int i, int j, double w0, double w1, double w3, double w4):
    return c_matrix_D_minus(i, j, w0, w1, w3, w4)

def flux_riemann(w_l, w_r):
    cdef double result[4]
    cdef double *_w_l
    cdef double *_w_r
    cdef int n
    numpy2c_double_inplace(w_l, &_w_l, &n)
    numpy2c_double_inplace(w_r, &_w_r, &n)
    c_flux_riemann(result, _w_l, _w_r)
    return c2numpy_double(&(result[0]), 4)

def flux_riemann_invert(w_l, w_r):
    cdef double result[4]
    cdef double *_w_l
    cdef double *_w_r
    cdef int n
    numpy2c_double_inplace(w_l, &_w_l, &n)
    numpy2c_double_inplace(w_r, &_w_r, &n)
    c_flux_riemann_invert(result, _w_l, _w_r)
    return c2numpy_double(&(result[0]), 4)

def numerical_flux(w_l, w_r, n):
    cdef double result[4]
    cdef double *_w_l
    cdef double *_w_r
    cdef double nx
    cdef double ny
    nx, ny = n
    cdef int _n
    numpy2c_double_inplace(w_l, &_w_l, &_n)
    numpy2c_double_inplace(w_r, &_w_r, &_n)
    c_numerical_flux(result, _w_l, _w_r, nx, ny)
    return c2numpy_double(&(result[0]), 4)
