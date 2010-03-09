#! /usr/bin/env python

# This example shows how to solve a first simple PDE:
#   - load the mesh,
#   - perform initial refinements
#   - create a H1 space over the mesh
#   - define weak formulation
#   - initialize matrix solver
#   - assemble and solve the matrix system
#   - visualize the solution
#
# PDE: Poisson equation -Laplace u = CONST_F with homogeneous (zero)
#      Dirichlet boundary conditions.
#
# You can change the constant right-hand side CONST_F, the
# initial polynomial degree P_INIT, and play with various initial
# mesh refinements at the beginning.

# Import modules
from hermes2d import Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space, \
        WeakForm, Solution, ScalarView, LinSystem, DummySolver
from hermes2d.forms import set_forms
from hermes2d.examples.c03 import set_bc
from hermes2d.examples import get_example_mesh

P_INIT = 5                # Uniform polynomial degree of mesh elements.

# Load the mesh file
mesh = Mesh()
mesh.load(get_example_mesh())

# Sample "manual" mesh refinement
#mesh.refine_element(0)

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
#sys.assemble()
#A = sys.get_matrix()
#b = sys.get_rhs()
#from scipy.sparse.linalg import cg
#x, res = cg(A, b)
#sln = Solution()
#sln.set_fe_solution(space, pss, x)

sln = Solution()
sys.assemble()
sys.solve_system(sln)

# Visualize the solution
view = ScalarView("Solution")
view.show(sln, lib="mayavi")

# Visualize the mesh
mview = MeshView("Hello world!", 100, 100, 500, 500)
mview.show(mesh, lib="mpl", method="orders", notebook=False)
