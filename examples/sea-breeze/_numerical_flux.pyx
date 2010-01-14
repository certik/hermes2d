from hermes2d._hermes2d cimport H1Space, c_H1Space, c_WeakForm, c_Solution, \
        WeakForm, Solution, c_Mesh, Mesh


cdef extern from "numerical_flux.h":
    double c_matrix_R "matrix_R"(int i, int j, double w0, double w1, double w3, double w4)
    double c_matrix_R_inv "matrix_R_inv"(int i, int j, double w0, double w1, double w3, double w4)

def matrix_R(int i, int j, double w0, double w1, double w3, double w4):
    return c_matrix_R(i, j, w0, w1, w3, w4)

def matrix_R_inv(int i, int j, double w0, double w1, double w3, double w4):
    return c_matrix_R_inv(i, j, w0, w1, w3, w4)
