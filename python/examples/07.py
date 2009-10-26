#! /usr/bin/env python

from hermes2d import Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space, \
        LinSystem, WeakForm, DummySolver, Solution, ScalarView, VonMisesFilter

from hermes2d.examples.c07 import set_bc, set_forms
from hermes2d.examples import get_sample_mesh

mesh = Mesh()
mesh.load(get_sample_mesh())
#mesh.refine_element(0)
#mesh.refine_all_elements()
#mesh.refine_towards_boundary(5, 3)
shapeset = H1Shapeset()
pss = PrecalcShapeset(shapeset)

# create an H1 space
xdisp = H1Space(mesh, shapeset)
ydisp = H1Space(mesh, shapeset)
xdisp.set_uniform_order(8)
ydisp.set_uniform_order(8)

set_bc(xdisp, ydisp)

ndofs = xdisp.assign_dofs(0)
ndofs += ydisp.assign_dofs(ndofs)

# initialize the discrete problem
wf = WeakForm(2)
set_forms(wf)

solver = DummySolver()
sys = LinSystem(wf, solver)
sys.set_spaces(xdisp, ydisp)
sys.set_pss(pss)

xsln = Solution()
ysln = Solution()
sys.assemble()
sys.solve_system(xsln, ysln)

view = ScalarView("Von Mises stress [Pa]", 50, 50, 1200, 600)
E = float(200e9)
nu = 0.3
stress = VonMisesFilter(xsln, ysln, E / (2*(1 + nu)),
        (E * nu) / ((1 + nu) * (1 - 2*nu)))
view.show(stress)

# view.wait()

mview = MeshView("Hello world!", 100, 100, 500, 500)
mview.show(mesh, lib="mpl", method="orders", notebook=False)
mview.wait()
