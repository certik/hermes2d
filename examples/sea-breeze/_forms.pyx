from hermes2d._hermes2d cimport H1Space, c_H1Space, c_WeakForm, c_Solution, \
        WeakForm, Solution, c_Mesh, Mesh


cdef extern from "forms.h":
    double c_TAU "TAU"
    double c_R "R"
    double c_c_v "c_v"
    void c_register_bc "register_bc"(c_H1Space s0, c_H1Space s1, c_H1Space s3,
            c_H1Space s4)
    void c_set_ic "set_ic"(c_Mesh mesh, c_Solution w0_prev,
            c_Solution w1_prev, c_Solution w3_prev, c_Solution w4_prev)
    void c_register_forms "register_forms"(c_WeakForm wf, c_Solution w0_prev,
            c_Solution w1_prev, c_Solution w3_prev, c_Solution w4_prev)

    double c_f_x "f_x"(int i, double w0, double w1, double w3, double w4)
    double c_f_z "f_z"(int i, double w0, double w1, double w3, double w4)
    double c_A_x "A_x"(int i, int j, double w0, double w1, double w3, double w4)
    double c_A_z "A_z"(int i, int j, double w0, double w1, double w3, double w4)

tau = c_TAU
R = c_R
c_v = c_c_v

def register_bc(H1Space s0, H1Space s1, H1Space s3, H1Space s4):
    c_register_bc((s0.thisptr)[0], (s1.thisptr)[0], (s3.thisptr)[0],
            (s4.thisptr)[0])

def set_ic(Mesh mesh, Solution w0_prev, Solution w1_prev, Solution
        w3_prev, Solution w4_prev):
    c_set_ic(mesh.thisptr[0],
            (<c_Solution*>w0_prev.thisptr)[0],
            (<c_Solution*>w1_prev.thisptr)[0],
            (<c_Solution*>w3_prev.thisptr)[0],
            (<c_Solution*>w4_prev.thisptr)[0]
            )

def register_forms(WeakForm wf, Solution w0_prev, Solution w1_prev, Solution
        w3_prev, Solution w4_prev):
    c_register_forms(wf.thisptr[0],
            (<c_Solution*>w0_prev.thisptr)[0],
            (<c_Solution*>w1_prev.thisptr)[0],
            (<c_Solution*>w3_prev.thisptr)[0],
            (<c_Solution*>w4_prev.thisptr)[0]
            )

def f_x(int i, double w0, double w1, double w3, double w4):
    return c_f_x(i, w0, w1, w3, w4)

def f_z(int i, double w0, double w1, double w3, double w4):
    return c_f_z(i, w0, w1, w3, w4)

def A_x(int i, int j, double w0, double w1, double w3, double w4):
    return c_A_x(i, j, w0, w1, w3, w4)

def A_z(int i, int j, double w0, double w1, double w3, double w4):
    return c_A_z(i, j, w0, w1, w3, w4)
