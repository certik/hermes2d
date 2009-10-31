#! /usr/bin/env python

from hermes2d import (Mesh, H1Shapeset, PrecalcShapeset, H1Space, Solution,
        WeakForm, DummySolver, LinSystem)
from _forms import register_bc, set_ic, register_forms

tau = 0.5

mesh = Mesh()
mesh.load("domain-quad.mesh")

mesh.refine_all_elements()
mesh.refine_all_elements()
mesh.refine_all_elements()

shapeset_h1 = H1Shapeset()
pss_h1 = PrecalcShapeset(shapeset_h1)
s0 = H1Space(mesh, shapeset_h1)
s1 = H1Space(mesh, shapeset_h1)
s3 = H1Space(mesh, shapeset_h1)
s4 = H1Space(mesh, shapeset_h1)

register_bc(s0, s1, s3, s4)

s0.set_uniform_order(1)
s1.set_uniform_order(1)
s3.set_uniform_order(1)
s4.set_uniform_order(1)

ndofs = 0
ndofs += s0.assign_dofs(ndofs)
ndofs += s1.assign_dofs(ndofs)
ndofs += s3.assign_dofs(ndofs)
ndofs += s4.assign_dofs(ndofs)

w0_prev = Solution()
w1_prev = Solution()
w3_prev = Solution()
w4_prev = Solution()
set_ic(mesh, w0_prev, w1_prev, w3_prev, w4_prev)

wf = WeakForm(4)
register_forms(wf, w0_prev, w1_prev, w3_prev, w4_prev)

solver = DummySolver()
sys = LinSystem(wf, solver)
sys.set_spaces(s0, s1, s3, s4)
sys.set_pss(pss_h1, pss_h1, pss_h1, pss_h1)


time = 0
for i in range(10):
    time += tau
    print "---- Time step %d, time = %g -----------------------------------" % \
            (i, time)

    ndofs = 0
    ndofs += s0.assign_dofs(ndofs)
    ndofs += s1.assign_dofs(ndofs)
    ndofs += s3.assign_dofs(ndofs)
    ndofs += s4.assign_dofs(ndofs)
    sys.assemble()
    w0_sln = Solution()
    w1_sln = Solution()
    w3_sln = Solution()
    w4_sln = Solution()
    sys.solve_system(w0_sln, w1_sln, w3_sln, w4_sln)
    #A = sys.get_matrix()
    #b = sys.get_rhs()
    w0_prev = w0_sln
    w1_prev = w1_sln
    w3_prev = w3_sln
    w4_prev = w4_sln
