#! /usr/bin/env python

# This example shows how to define Neumann boundary conditions. In addition,
# you will see how a Filter is used to visualize gradient of the solution
#
# PDE: Poisson equation -Laplace u = f, where f = CONST_F
#
# BC: u = 0 on Gamma_4 (edges meeting at the re-entrant corner)
#     du/dn = CONST_GAMMA_1 on Gamma_1 (bottom edge)
#     du/dn = CONST_GAMMA_2 on Gamma_2 (top edge, circular arc, and right-most edge)
#     du/dn = CONST_GAMMA_3 on Gamma_3 (left-most edge)
#
# You can play with the parameters below. For most choices of the four constants,
# the solution has a singular (infinite) gradient at the re-entrant corner.
# Therefore we visualize not only the solution but also its gradient.

# Import modules
from hermes2d import Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space, \
        LinSystem, Solution, ScalarView, WeakForm, DummySolver

from hermes2d.examples.c05 import set_bc, set_forms
from hermes2d.examples.c05 import set_forms as set_forms_surf
from hermes2d.forms import set_forms
from hermes2d.examples import get_example_mesh

CONST_F = -4.0                       # right-hand side
P_INIT = 4                           # initial polynomial degree in all elements
CORNER_REF_LEVEL = 12                # number of mesh refinements towards the re-entrant corner

# Load the mesh file
mesh = Mesh()
mesh.load(get_example_mesh())
#mesh.refine_towards_vertex(3, CORNER_REF_LEVEL)

# Initialize the shapeset and the cache
shapeset = H1Shapeset()
pss = PrecalcShapeset(shapeset)

# Create an H1 space
space = H1Space(mesh, shapeset)
space.set_uniform_order(P_INIT)
set_bc(space)
space.assign_dofs()

# Initialize the discrete problem
wf = WeakForm(1)
set_forms(wf, CONST_F)
set_forms_surf(wf)

# Assemble the stiffness matrix and solve the system
sln = Solution()
solver = DummySolver()
sys = LinSystem(wf, solver)
sys.set_spaces(space)
sys.set_pss(pss)
sys.assemble()
sys.solve_system(sln)

# Visualize the approximation
view = ScalarView("Solution")
view.show(sln, lib="mayavi")

# Visualize the mesh
mview = MeshView("Hello world!", 100, 100, 500, 500)
mview.show(mesh, lib="mpl", method="orders", notebook=False)
