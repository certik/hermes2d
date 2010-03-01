#! /usr/bin/env python

# This example shows how to solve a time-dependent PDE discretized
# in time via the implicit Euler method. The St. Vitus Cathedral
# in Prague (http://en.wikipedia.org/wiki/St._Vitus_Cathedral)
# responds to changes in the surrounding air temperature
# during one 24-hour cycle. You will also learn how to use the
# solution calculated in the previous time step.
#
# PDE: non-stationary heat transfer equation
# HEATCAP * RHO * dT/dt - LAMBDA * Laplace T = 0
#
# Domain: St. Vitus cathedral (cathedral.mesh)
#
# IC:  T = T_INIT
# BC:  T = T_INIT on the bottom edge ... Dirichlet
#      dT/dn = ALPHA*(t_exterior(time) - T) ... Newton, time-dependent
#
# Time-stepping: implicit Euler

# Import modules
from hermes2d import Mesh, MeshView, VectorView, OrderView, H1Shapeset, PrecalcShapeset, H1Space, \
WeakForm, Solution, ScalarView, LinSystem, DummySolver, RefSystem, H1OrthoHP

from hermes2d.examples.c09 import set_bc, temp_ext, set_forms, update_time
from hermes2d.examples import get_cathedral_mesh

# The following parameters can be played with:
P_INIT = 1            # polynomial degree of elements
INIT_REF_NUM = 4      # number of initial uniform refinements
TAU = 300.0           # time step in seconds

# Problem constants
T_INIT = 10           # temperature of the ground (also initial temperature)
FINAL_TIME = 86400    # length of time interval (24 hours) in seconds

# Global variable
TIME = 0;

# Load the mesh
mesh = Mesh()
mesh.load(get_cathedral_mesh())

for i in range(INIT_REF_NUM):
    mesh.refine_all_elements()
mesh.refine_towards_boundary(2, 5)

# Set up shapeset
shapeset = H1Shapeset()
pss = PrecalcShapeset(shapeset)

# Set up spaces
space = H1Space(mesh, shapeset)
set_bc(space)
space.set_uniform_order(P_INIT)

# Enumerate basis functions
space.assign_dofs()

# Set initial condition
tsln = Solution()
tsln.set_const(mesh, T_INIT)

# Weak formulation
wf = WeakForm(1)
set_forms(wf, tsln)

# Matrix solver
solver = DummySolver()

# Linear system
ls = LinSystem(wf, solver)
ls.set_spaces(space)
ls.set_pss(pss)

# Visualisation
sview = ScalarView("Temperature", 0, 0, 450, 600)
#title = "Time %s, exterior temperature %s" % (TIME, temp_ext(TIME))
#Tview.set_min_max_range(0,20);
#Tview.set_title(title);
#Tview.fix_scale_width(3);

# Time stepping
nsteps = int(FINAL_TIME/TAU + 0.5)
rhsonly = False;

for n in range(1,nsteps+1):
    print ("\n---- Time %s, time step %s, ext_temp %s ----------" % (TIME, n, temp_ext(TIME)) )

    # Assemble and solve
    ls.assemble()
    rhsonly = True
    ls.solve_system(tsln, lib="scipy")

    # Shifting the time variable
    TIME += TAU
    update_time(TIME)

    # Visualization of solution
    title = "Time %s, exterior temperature %s" % (TIME, temp_ext(TIME))
    sview.show(tsln, lib="mayavi")
