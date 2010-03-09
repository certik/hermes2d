#! /usr/bin/env python

# This example explains how to use Newton boundary conditions. Again,
# a Filter is used to visualize the solution gradient.
#
# PDE: Laplace equation -Laplace u = 0 (this equation describes, among
# many other things, also stationary heat transfer in a homogeneous linear
# material).
#
# BC: u = T1 ... fixed temperature on Gamma_3 (Dirichlet)
#     du/dn = 0 ... insulated wall on Gamma_2 and Gamma_4 (Neumann)
#     du/dn = H*(u - T0) ... heat flux on Gamma_1 (Newton)
#
# Note that the last BC can be written in the form  du/dn - H*u = -h*T0.

# Import modules
from hermes2d import Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space, \
        LinSystem, WeakForm, DummySolver, Solution, ScalarView

from hermes2d.examples.c06 import set_bc, set_forms
from hermes2d.examples import get_example_mesh

# You can play with the parameters below:
T1 = 30.0                # prescribed temperature on Gamma_3
T0 = 20.0                # outer temperature on Gamma_1
H  = 0.05                # heat flux on Gamma_1
P_INIT = 6               # uniform polynomial degree in the mesh
UNIFORM_REF_LEVEL = 2    # number of initial uniform mesh refinements
CORNER_REF_LEVEL = 12    # number of mesh refinements towards the re-entrant corner

# Load the mesh file
mesh = Mesh()
mesh.load(get_example_mesh())
#for i in range(UNIFORM_REF_LEVEL):
#    mesh.refine_all_elemenrs()
#mesh.refine_towards_vertex(3, CORNER_REF_LEVEL)

# Initialize the shapeset and the cache
shapeset = H1Shapeset()
pss = PrecalcShapeset(shapeset)

# Create an H1 space
space = H1Space(mesh, shapeset)
space.set_uniform_order(P_INIT)
set_bc(space)
space.assign_dofs()

# Initialize the weak formulation
wf = WeakForm(1)
set_forms(wf)

# Initialize the linear system and solver
solver = DummySolver()
sys = LinSystem(wf, solver)
sys.set_spaces(space)
sys.set_pss(pss)

# Assemble the stiffness matrix and solve the system
sln = Solution()
sys.assemble()
sys.solve_system(sln)

# Visualize the solution
view = ScalarView("Solution")
view.show(sln, lib="mayavi")

# Visualize the mesh
mview = MeshView("Hello world!", 100, 100, 500, 500)
mview.show(mesh, lib="mpl", method="orders", notebook=False)
