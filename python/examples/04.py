#! /usr/bin/env python

from hermes2d import (Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space,
        LinSystem, Solution, ScalarView, WeakForm, DummySolver)

from hermes2d.examples.c04 import set_bc
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
wf = WeakForm()
set_forms(wf, -4)

solver = DummySolver()
sys = LinSystem(wf, solver)
sys.set_spaces(space)
sys.set_pss(pss)

# assemble the stiffness matrix and solve the system
sys.assemble()
sln = Solution()
sys.solve_system(sln)

view = ScalarView("Solution")
view.show(sln, lib="mayavi")
# view.wait()

mview = MeshView("Hello world!", 100, 100, 500, 500)
mview.show(mesh, lib="mpl", method="orders", notebook=False)
mview.wait()
