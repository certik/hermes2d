#! /usr/bin/env python

# This example explains how to create two spaces over a mesh and use them
# to solve a simple problem of linear elasticity. At the end, VonMises
# filter is used to visualize the stress.
#
# PDE: Lame equations of linear elasticity
#
# BC: du_1/dn = f_0 on Gamma_3 and du_1/dn = 0 on Gamma_2, Gamma_4, Gamma_5
#     du_2/dn = f_1 on Gamma_3 and du_2/dn = 0 on Gamma_2, Gamma_4, Gamma_5
#     u_1 = 0 and u_2 = 0 on Gamma_1

# Import modules
from hermes2d import Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space, \
        LinSystem, WeakForm, DummySolver, Solution, ScalarView, VonMisesFilter

from hermes2d.examples.c08 import set_bc, set_forms
from hermes2d.examples import get_sample_mesh

# The following parameter can be changed:
P_INIT = 8

# Load the mesh file
mesh = Mesh()
mesh.load(get_sample_mesh())

# Initialize the shapeset and the cache
shapeset = H1Shapeset()
pss = PrecalcShapeset(shapeset)

# Create the x displacement space
xdisp = H1Space(mesh, shapeset)
xdisp.set_uniform_order(P_INIT)

# Create the y displacement space
ydisp = H1Space(mesh, shapeset)
ydisp.set_uniform_order(P_INIT)

set_bc(xdisp, ydisp)
ndofs = xdisp.assign_dofs(0)
ndofs += ydisp.assign_dofs(ndofs)

# Initialize the weak formulation
wf = WeakForm(2)
set_forms(wf)

# Initialize the linear system and solver
solver = DummySolver()
sys = LinSystem(wf, solver)
sys.set_spaces(xdisp, ydisp)
sys.set_pss(pss)

# Assemble the stiffness matrix and solve the system
xsln = Solution()
ysln = Solution()
sys.assemble()
sys.solve_system(xsln, ysln, lib="scipy")

# Visualize the solution
view = ScalarView("Von Mises stress [Pa]", 50, 50, 1200, 600)
E = float(200e9)
nu = 0.3
l = (E * nu) / ((1 + nu) * (1 - 2*nu))
mu = E / (2*(1 + nu))
stress = VonMisesFilter(xsln, ysln, mu, l)
view.show(stress, lib="mayavi")

# Visualize the mesh
mview = MeshView("Hello world!", 100, 100, 500, 500)
mview.show(mesh, lib="mpl", method="orders", notebook=False)
