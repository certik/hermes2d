#! /usr/bin/env python

# This example explains how to use the multimesh adaptive hp-FEM,
# where different physical fields (or solution components) can be
# approximated using different meshes and equipped with mutually
# independent adaptivity mechanisms. Here we consider linear elasticity
# and will approximate each displacement components using an individual
# mesh.
#
# PDE: Lame equations of linear elasticity, treated as a coupled system
#      of two PDEs
#
# BC: u_1 = u_2 = 0 on Gamma_1
#     du_2/dn = f on Gamma_2
#     du_1/dn = du_2/dn = 0 elsewhere

# Import modules
from hermes2d import Mesh, MeshView, VectorView, OrderView, H1Shapeset, PrecalcShapeset, H1Space, \
        WeakForm, Solution, ScalarView, LinSystem, DummySolver, RefSystem, \
        H1Adapt, H1ProjBasedSelector, CandList, \
        VonMisesFilter

from hermes2d.examples.c11 import set_bc, set_wf_forms, set_hp_forms
from hermes2d.examples import get_bracket_mesh

# The following parameters can be changed: In particular, compare hp- and
# h-adaptivity via the ADAPT_TYPE option, and compare the multi-mesh vs. single-mesh
# using the MULTI parameter.
P_INIT = 1               # Initial polynomial degree of all mesh elements.
MULTI = True             # MULTI = true  ... use multi-mesh,
                            # MULTI = false ... use single-mesh.
                            # Note: In the single mesh option, the meshes are
                            # forced to be geometrically the same but the
                            # polynomial degrees can still vary.
SAME_ORDERS = True       # SAME_ORDERS = true ... when single-mesh is used,
                            # this forces the meshes for all components to be
                            # identical, including the polynomial degrees of
                            # corresponding elements. When multi-mesh is used,
                            # this parameter is ignored.
THRESHOLD = 0.3          # This is a quantitative parameter of the adapt(...) function and
                                 # it has different meanings for various adaptive strategies (see below).
STRATEGY = 1             # Adaptive strategy:
                            # STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                            #   error is processed. If more elements have similar errors, refine
                            #   all to keep the mesh symmetric.
                            # STRATEGY = 1 ... refine all elements whose error is larger
                            #   than THRESHOLD times maximum element error.
                            # STRATEGY = 2 ... refine all elements whose error is larger
                            #   than THRESHOLD.
                            # More adaptive strategies can be created in adapt_ortho_h1.cpp.

CAND_TYPE = CandList.HP_ANISO  # Predefined list of element refinement candidates.
                        # Possible values are are attributes of the class CandList:
                        # P_ISO, P_ANISO, H_ISO, H_ANISO, HP_ISO, HP_ANISO_H, HP_ANISO_P, HP_ANISO
                        # See the Sphinx tutorial (http://hpfem.org/hermes2d/doc/src/tutorial-2.html#adaptive-h-fem-and-hp-fem) for details.

MESH_REGULARITY = -1     # Maximum allowed level of hanging nodes:
                            # MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                            # MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                            # MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                            # Note that regular meshes are not supported, this is due to
                            # their notoriously bad performance.
MAX_ORDER = 10           # Maximum allowed element degree
ERR_STOP = 0.5           # Stopping criterion for adaptivity (rel. error tolerance between the
                            # fine mesh and coarse mesh solution in percent).
NDOF_STOP = 40000        # Adaptivity process stops when the number of degrees of freedom grows over
                            # this limit. This is mainly to prevent h-adaptivity to go on forever.

# Problem constants
E  = 200e9               # Young modulus for steel: 200 GPa
nu = 0.3                 # Poisson ratio
lamda = (E * nu) / ((1 + nu) * (1 - 2*nu))
mu = E / (2*(1 + nu))

# Load the mesh
xmesh = Mesh()
ymesh = Mesh()
xmesh.load(get_bracket_mesh())

# Create initial mesh for the vertical displacement component,
# identical to the mesh for the horizontal displacement
# (bracket.mesh becomes a master mesh)
ymesh.copy(xmesh)

# Initialize the shapeset and the cache
shapeset = H1Shapeset()
xpss = PrecalcShapeset(shapeset)
ypss = PrecalcShapeset(shapeset)

# Create the x displacement space
xdisp = H1Space(xmesh, shapeset)
set_bc(xdisp)
xdisp.set_uniform_order(P_INIT)

# Create the x displacement space
ydisp = H1Space(ymesh, shapeset)
set_bc(ydisp)
ydisp.set_uniform_order(P_INIT)

# Enumerate basis functions
ndofs = xdisp.assign_dofs()
ydisp.assign_dofs(ndofs)

# Initialize the weak formulation
wf = WeakForm(2)
set_wf_forms(wf)

# Visualization of solution and meshes
xoview = OrderView("X polynomial orders", 0, 0, 500, 500)
yoview = OrderView("Y polynomial orders", 510, 0, 500, 500)
sview = ScalarView("Von Mises stress [Pa]", 1020, 0, 500, 500)

# Matrix solver
solver = DummySolver()

# adaptivity loop
it = 1
done = False
cpu = 0.0

x_sln_coarse = Solution()
y_sln_coarse = Solution()

x_sln_fine = Solution()
y_sln_fine = Solution()

selector = H1ProjBasedSelector(CAND_TYPE, 1.0, MAX_ORDER, shapeset)

while(not done):

    print ("\n---- Adaptivity step %d ---------------------------------------------\n" % it)
    it += 1

    # Calculating the number of degrees of freedom
    ndofs = xdisp.assign_dofs()
    ndofs += ydisp.assign_dofs(ndofs)

    print("xdof=%d, ydof=%d\n" % (xdisp.get_num_dofs(), ydisp.get_num_dofs()) )

    # Solve the coarse mesh problem
    ls = LinSystem(wf, solver)
    ls.set_spaces(xdisp, ydisp)
    ls.set_pss(xpss, ypss)
    ls.assemble()
    ls.solve_system(x_sln_coarse, y_sln_coarse, lib="scipy")

    # View the solution -- this can be slow; for illustration only
    stress_coarse = VonMisesFilter(x_sln_coarse, y_sln_coarse, mu, lamda)
    #sview.set_min_max_range(0, 3e4)
    sview.show(stress_coarse)
    #xoview.show(xdisp)
    #yoview.show(ydisp)
    xmesh.plot(space=xdisp)
    ymesh.plot(space=ydisp)

    # Solve the fine mesh problem
    rs = RefSystem(ls)
    rs.assemble()
    rs.solve_system(x_sln_fine, y_sln_fine, lib="scipy")

    # Calculate element errors and total error estimate
    hp = H1Adapt([xdisp, ydisp])
    hp.set_solutions([x_sln_coarse, y_sln_coarse], [x_sln_fine, y_sln_fine]);
    set_hp_forms(hp)
    err_est = hp.calc_error() * 100

    print("Error estimate: %s" % err_est)

# If err_est too large, adapt the mesh
    if err_est < ERR_STOP:
        done = True
    else:
        hp.adapt(selector, THRESHOLD, STRATEGY, MESH_REGULARITY, SAME_ORDERS)
        ndofs = xdisp.assign_dofs()
        ndofs += ydisp.assign_dofs(ndofs)

        if ndofs >= NDOF_STOP:
            done = True


# Show the fine solution - this is the final result
stress_fine = VonMisesFilter(x_sln_fine, y_sln_fine, mu, lamda)
sview.show(stress_fine)
