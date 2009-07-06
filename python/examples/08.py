#! /usr/bin/env python

from hermes2d import Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
        LinSystem, WeakForm, DummySolver, Solution, ScalarView, VectorView

from hermes2d.examples.c08 import set_bc, set_forms
from hermes2d.examples import get_cylinder_mesh

mesh = Mesh()
mesh.load(get_cylinder_mesh())
#mesh.refine_element(0)
#mesh.refine_all_elements()
mesh.refine_towards_boundary(5, 3)
shapeset = H1Shapeset()
pss = PrecalcShapeset(shapeset)

# create an H1 space
xvel = H1Space(mesh, shapeset)
yvel = H1Space(mesh, shapeset)
press = H1Space(mesh, shapeset)
xvel.set_uniform_order(2)
yvel.set_uniform_order(2)
press.set_uniform_order(1)

set_bc(xvel, yvel, press)

ndofs = 0
ndofs += xvel.assign_dofs(ndofs)
ndofs += yvel.assign_dofs(ndofs)
ndofs += press.assign_dofs(ndofs)

xprev = Solution()
yprev = Solution()

xprev.set_zero(mesh)
yprev.set_zero(mesh)

# initialize the discrete problem
wf = WeakForm(3)
set_forms(wf, xprev, yprev)

# visualize the solution
vview = VectorView("velocity [m/s]", 0, 0, 1200, 350)
pview = ScalarView("pressure [Pa]", 0, 500, 1200, 350)
vview.set_min_max_range(0, 1.9)
vview.show_scale(False)
pview.show_scale(False)
pview.show_mesh(False)

solver = DummySolver()
sys = LinSystem(wf, solver)
sys.set_spaces(xvel, yvel, press)
sys.set_pss(pss)
#dp.set_external_fns(xprev, yprev)

EPS_LOW = 0.0014

for i in range(1000):
    print "*** Iteration %d ***" % i
    psln = Solution()
    sys.assemble()
    sys.solve_system(xprev, yprev, psln)
    vview.show(xprev, yprev, 2*EPS_LOW)
    pview.show(psln)

vview.wait()
