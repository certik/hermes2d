#! /usr/bin/env python

# This example illustrates how to use non-homogeneous(nonzero)
# Dirichlet boundary conditions.
#
# PDE: Poisson equation -Laplace u = CONST_F, where CONST_F is
# a constant right-hand side. It is not difficult to see that
# the function u(x,y) = (-CONST_F/4)*(x^2 + y^2) satisfies the
# above PDE. Since also the Dirichlet boundary conditions
# are chosen to match u(x,y), this function is the exact
# solution.
# 
# Note that since the exact solution is a quadratic polynomial,
# Hermes will compute it exactly if all mesh elements are quadratic
# or higher (then the exact solution lies in the finite element space).
# If some elements in the mesh are linear, Hermes will only find
# an approximation.

# Import modules
from hermes2d import (Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space,
        LinSystem, Solution, ScalarView, WeakForm, DummySolver)

from hermes2d.examples.c04 import set_bc
from hermes2d.examples import get_example_mesh
from hermes2d.forms import set_forms

# Below you can play with the parameters CONST_F, P_INIT, and UNIFORM_REF_LEVEL.
CONST_F = -4.0           # constant right-hand side
P_INIT = 2               # initial polynomial degree in all elements
UNIFORM_REF_LEVEL = 3    # number of initial uniform mesh refinements

# Load the mesh file
mesh = Mesh()
mesh.load(get_example_mesh())
#for i in range(UNIFORM_REF_LEVEL):
#    mesh.refine_all_elemenrs()

# Initialize the shapeset and the cache
shapeset = H1Shapeset()
pss = PrecalcShapeset(shapeset)

# Create an H1 space
space = H1Space(mesh, shapeset)
space.set_uniform_order(P_INIT)
set_bc(space)
space.assign_dofs()

# Initialize the weak formulation
wf = WeakForm()
set_forms(wf, CONST_F)

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
