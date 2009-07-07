#! /usr/bin/env python

from hermes2d import Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
        LinSystem, Solution, ScalarView, WeakForm, DummySolver

from hermes2d.examples.c05 import set_bc, set_forms
from hermes2d.examples.c05 import set_forms as set_forms_surf
from hermes2d.forms import set_forms
from hermes2d.examples import get_example_mesh

mesh = Mesh()
mesh.load(get_example_mesh())
#mesh.refine_element(0)
#mesh.refine_all_elements()
mesh.refine_towards_boundary(5, 3)
shapeset = H1Shapeset()
pss = PrecalcShapeset(shapeset)

# create an H1 space
space = H1Space(mesh, shapeset)
space.set_uniform_order(5)

set_bc(space)

space.assign_dofs()

xprev = Solution()
yprev = Solution()

# initialize the discrete problem
wf = WeakForm(1)
set_forms(wf, -1)
set_forms_surf(wf)

sln = Solution()
solver = DummySolver()
sys = LinSystem(wf, solver)
sys.set_spaces(space)
sys.set_pss(pss)
sys.assemble()
sys.solve_system(sln)

view = ScalarView("Solution")
view.show(sln)
view.wait()
