from hermes2d._hermes2d cimport H1Space, c_H1Space, c_WeakForm, c_Solution, \
        WeakForm, Solution, c_Mesh, Mesh
from hermes2d._hermes2d cimport array_double_c2numpy, \
        array_double_numpy2c_inplace


cdef extern from "numerical_flux.h":
    double c_matrix_R "matrix_R"(int i, int j, double w0, double w1, double w3, double w4)
    double c_matrix_R_inv "matrix_R_inv"(int i, int j, double w0, double w1, double w3, double w4)
    double c_matrix_D_minus "matrix_D_minus"(int i, int j, double w0, double w1, double w3, double w4)
    void c_flux_riemann "flux_riemann"(double result[4], double w_l[4], double w_r[4])
    double c_R "R"
    double c_c_v "c_v"

R = c_R
c_v = c_c_v

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
    array_double_numpy2c_inplace(w_l, &_w_l, &n)
    array_double_numpy2c_inplace(w_r, &_w_r, &n)
    c_flux_riemann(result, _w_l, _w_r)
    return array_double_c2numpy(&(result[0]), 4)
