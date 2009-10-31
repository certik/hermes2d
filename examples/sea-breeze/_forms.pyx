from hermes2d._hermes2d cimport H1Space, c_H1Space, c_WeakForm, c_Solution, \
        WeakForm, Solution, c_Mesh, Mesh

cdef extern from "forms.h":
    void c_register_bc "register_bc"(c_H1Space s0, c_H1Space s1, c_H1Space s3,
            c_H1Space s4)
    void c_set_ic "set_ic"(c_Mesh mesh, c_Solution w0_prev,
            c_Solution w1_prev, c_Solution w3_prev, c_Solution w4_prev)
    void c_register_forms "register_forms"(c_WeakForm wf, c_Solution w0_prev,
            c_Solution w1_prev, c_Solution w3_prev, c_Solution w4_prev)

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
