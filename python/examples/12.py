#! /usr/bin/env python

#  This example uses automatic adaptivity to solve a general second-order linear
#  equation with non-constant coefficients.
#
#  PDE: -d/dx(a_11(x,y)du/dx) - d/dx(a_12(x,y)du/dy) - d/dy(a_21(x,y)du/dx) - d/dy(a_22(x,y)du/dy)
#       + a_1(x,y)du/dx + a_21(x,y)du/dy + a_0(x,y)u = rhs(x,y)
#
#  Domain: arbitrary
#
#  BC:  Dirichlet for boundary marker 1: u = g_D(x,y)
#       Natural for any other boundary marker:   (a_11(x,y)*nu_1 + a_21(x,y)*nu_2) * dudx
#                                              + (a_12(x,y)*nu_1 + s_22(x,y)*nu_2) * dudy = g_N(x,y)

# Import modules
from hermes2d import Mesh, MeshView, VectorView, OrderView, H1Shapeset, PrecalcShapeset, H1Space, \
WeakForm, Solution, ScalarView, LinSystem, DummySolver, RefSystem, H1OrthoHP

from hermes2d.examples.c12 import set_bc, set_forms
from hermes2d.examples import get_12_mesh

#  The following parameters can be changed:
P_INIT = 1              # Initial polynomial degree of all mesh elements.
THRESHOLD = 0.6         # This is a quantitative parameter of the adapt(...) function and
                        # it has different meanings for various adaptive strategies (see below).
STRATEGY = 0            # Adaptive strategy:
                            # STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                            #   error is processed. If more elements have similar errors, refine
                            #   all to keep the mesh symmetric.
                            # STRATEGY = 1 ... refine all elements whose error is larger
                            #   than THRESHOLD times maximum element error.
                            # STRATEGY = 2 ... refine all elements whose error is larger
                            #   than THRESHOLD.
                            # More adaptive strategies can be created in adapt_ortho_h1.cpp.
ADAPT_TYPE = 0          # Type of automatic adaptivity:
                            # ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                            # ADAPT_TYPE = 1 ... adaptive h-FEM,
                            # ADAPT_TYPE = 2 ... adaptive p-FEM.
ISO_ONLY = False        # Isotropic refinement flag (concerns quadrilateral elements only).
                            # ISO_ONLY = false ... anisotropic refinement of quad elements
                            # is allowed (default),
                            # ISO_ONLY = true ... only isotropic refinements of quad elements
                            # are allowed.
MESH_REGULARITY = -1    # Maximum allowed level of hanging nodes:
                            # MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                            # MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                            # MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                            # Note that regular meshes are not supported, this is due to
                            # their notoriously bad performance.
ERR_STOP = 0.01         # Stopping criterion for adaptivity (rel. error tolerance between the
                            # fine mesh and coarse mesh solution in percent).
NDOF_STOP = 40000       # Adaptivity process stops when the number of degrees of freedom grows
                            # over this limit. This is to prevent h-adaptivity to go on forever.

# Load the mesh
mesh = Mesh()
mesh.load(get_12_mesh())
mesh.refine_all_elements()

# Initialize the shapeset and the cache
shapeset = H1Shapeset()
pss = PrecalcShapeset(shapeset)

# Create finite element space
space = H1Space(mesh, shapeset)
set_bc(space)
space.set_uniform_order(P_INIT)

# Enumerate basis functions
space.assign_dofs()

# Initialize the weak formulation
wf = WeakForm(1)
set_forms(wf)

# Visualize solution and mesh
sview = ScalarView("Coarse solution", 0, 0, 600, 1000)
oview = OrderView("Polynomial orders", 1220, 0, 600, 1000)
mview = MeshView("Example 12", 100, 100, 500, 500)

# Matrix solver
solver = DummySolver()

# Adaptivity loop
it = 0
ndofs = 0
done = False
sln_coarse = Solution()
sln_fine = Solution()

while (not done):
    print("\n---- Adaptivity step %d ---------------------------------------------\n" % (it+1))
    it += 1

    # Solve the coarse mesh problem
    ls = LinSystem(wf, solver)
    ls.set_spaces(space)
    ls.set_pss(pss)

    ls.assemble()
    ls.solve_system(sln_coarse)

    # View the solution and mesh
    sview.show(sln_coarse, lib='mayavi');
    mview.show(mesh, space=space, lib="mpl", method="orders", notebook=False)

    # Solve the fine mesh problem
    rs = RefSystem(ls)
    rs.assemble()
    rs.solve_system(sln_fine)

    # Calculate element errors and total error estimate
    hp = H1OrthoHP(space);
    err_est = hp.calc_error(sln_coarse, sln_fine) * 100
    print("Error estimate: %d" % err_est)

    # If err_est too large, adapt the mesh
    if (err_est < ERR_STOP):
        done = True
    else:
        hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE)#, ISO_ONLY, MESH_REGULARITY)
        ndofs = space.assign_dofs()

        if (ndofs >= NDOF_STOP):
            done = True

sview.show(sln_fine, lib='mayavi')
