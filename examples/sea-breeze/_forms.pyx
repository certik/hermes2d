from hermes2d._hermes2d cimport H1Space, c_H1Space

cdef extern from "forms.h":
    void c_register_bc "register_bc"(c_H1Space s0, c_H1Space s1, c_H1Space s3,
            c_H1Space s4)

def register_bc(H1Space s0, H1Space s1, H1Space s3, H1Space s4):
    c_register_bc((s0.thisptr)[0], (s1.thisptr)[0], (s3.thisptr)[0],
            (s4.thisptr)[0])
