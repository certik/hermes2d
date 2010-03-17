#! /usr/bin/env python

from hermes2d import (Mesh, MeshView, OrderView, H1Shapeset, L2Shapeset,
        PrecalcShapeset, H1Space, LinSystem, WeakForm, DummySolver, Solution,
        ScalarView, VonMisesFilter, OrderView)

from hermes2d.examples.c18 import set_bc, set_bc_x, set_bc_y, set_forms
from hermes2d.examples import get_18_mesh

# The time-dependent laminar incompressible Navier-Stokes equations are
# discretized in time via the implicit Euler method. If NEWTON == true,
# the Newton's method is used to solve the nonlinear problem at each time 
# step. If NEWTON == false, the convective term is only linearized using the 
# velocities from the previous time step. Obviously the latter approach is wrong, 
# but people do this frequently because it is faster and simpler to implement. 
# Therefore we include this case for comparison purposes. We also show how 
# to use discontinuous ($L^2$) elements for pressure and thus make the 
# velocity discreetely divergence free. Comparison to approximating the 
# pressure with the standard (continuous) Taylor-Hood elements is enabled.  
# The Reynolds number Re = 200 which is embarrassingly low. You 
# can increase it but then you will need to make the mesh finer, and the 
# computation will take more time. 
#
# PDE: incompressible Navier-Stokes equations in the form
# \partial v / \partial t - \Delta v / Re + (v \cdot \nabla) v + \nabla p = 0,
# div v = 0
#
# BC: u_1 is a time-dependent constant and u_2 = 0 on Gamma_4 (inlet)
#     u_1 = u_2 = 0 on Gamma_1 (bottom), Gamma_3 (top) and Gamma_5 (obstacle)
#     "do nothing" on Gamma_2 (outlet)
#
# Geometry: Rectangular channel containing an off-axis circular obstacle. The
#           radius and position of the circle, as well as other geometry 
#           parameters can be changed in the mesh file "domain.mesh". 
#
#
# The following parameters can be changed:
#

PRESSURE_IN_L2 = False
#NEWTON = True
NEWTON = False

P_INIT_VEL = 2
P_INIT_PRESSURE = 1

RE = 200.0
VEL_INLET = 1.0
STARTUP_TIME = 1.0

TAU = 0.1
T_FINAL = 30000.0
NEWTON_TOL = 1e-3
NEWTON_MAX_ITER = 10
H = 5

# Boundary markers in the mesh file
marker_bottom = 1
marker_right  = 2
marker_top = 3
marker_left = 4
marker_obstacle = 5

# Current time (defined as global since needed in weak forms)
TIME = 0

# Load the mesh file
mesh = Mesh()
mesh.load(get_18_mesh())

# A-priori mesh refinements
# mesh.refine_all_elements()
# mesh.refine_towards_boundary(5, 4, false)
# mesh.refine_towards_boundary(1, 4)
# mesh.refine_towards_boundary(3, 4)

# Initialize shapesets and the cache
h1_shapeset = H1Shapeset()
h1_pss = PrecalcShapeset(h1_shapeset)

if PRESSURE_IN_L2:
  l2_shapeset = L2Shapeset()
  l2_pss = PrecalcShapeset(l2_shapeset)

# Spaces for velocities and pressure
xvel_space = H1Space(mesh, h1_shapeset)
yvel_space = H1Space(mesh, h1_shapeset)

if PRESSURE_IN_L2:
    p_space = L2Space(mesh, l2_shapeset)
else:
    p_space = H1Space(mesh, h1_shapeset)

set_bc_x(xvel_space)
set_bc_y(yvel_space)
set_bc(p_space)

# Set velocity and pressure polynomial degrees
xvel_space.set_uniform_order(P_INIT_VEL)
yvel_space.set_uniform_order(P_INIT_VEL)
p_space.set_uniform_order(P_INIT_PRESSURE)

# Assign degrees of freedom
ndofs = 0
ndofs += xvel_space.assign_dofs(ndofs)
ndofs += yvel_space.assign_dofs(ndofs)
ndofs += p_space.assign_dofs(ndofs)

# Solutions for the Newton's iteration and time stepping
xvel_prev_time = Solution()
yvel_prev_time = Solution()
xvel_prev_newton = Solution()
p_prev = Solution()
yvel_prev_newton = Solution()

xvel_prev_time.set_zero(mesh)
yvel_prev_time.set_zero(mesh)
xvel_prev_newton.set_zero(mesh)
yvel_prev_newton.set_zero(mesh)
p_prev.set_zero(mesh)

# Set up weak formulation
wf = WeakForm(3)

if (NEWTON):
    set_form_n(wf)
else:
    set_forms(wf)

# Visualization
#vview = VectorView("velocity [m/s]", 0, 0, 1500, 470)
pview = ScalarView("pressure [Pa]", 0, 530, 1500, 470)
#vview.set_min_max_range(0, 1.6)
#vview.fix_scale_width(80)
#pview.set_min_max_range(-0.9, 1.0)
#pview.fix_scale_width(80)
pview.show_mesh(True)

# Matrix solver
umfpack = DummySolver()

# Linear system
ls = LinSystem(wf, umfpack)

# Nonlinear system
#nls = NonlinSystem(wf, umfpack)

if (NEWTON):
    # Set up the nonlinear system
    nls.set_spaces(xvel_space, yvel_space, p_space)
    if PRESSURE_IN_L2:
        nls.set_pss(h1_pss, h1_pss, l2_pss)
    else:
        nls.set_pss(h1_pss)
else:
    # Set up the linear system
    ls.set_spaces(xvel_space, yvel_space, p_space)
    if PRESSURE_IN_L2:
        ls.set_pss(h1_pss, h1_pss, l2_pss)
    else:
        ls.set_pss(h1_pss)

# Time-stepping loop
title = [0]*100
num_time_steps = T_FINAL / TAU

for i in range(1, num_time_steps):
    TIME += TAU

    print("\n---- Time step %d, time = %g:\n"% (i, TIME))
    
    # This is needed to update the time-dependent boundary conditions
    ndofs = 0
    ndofs += xvel_space.assign_dofs(ndofs)
    ndofs += yvel_space.assign_dofs(ndofs)
    ndofs += p_space.assign_dofs(ndofs)

    if (NEWTON):
        # Newton's method
        if (not nls.solve_newton_3(xvel_prev_newton, yvel_prev_newton, p_prev, NEWTON_TOL, NEWTON_MAX_ITER)):
            error("Newton's method did not converge.")
        # Show the solution at the end of time step
        sprintf(title, "Velocity, time %g", TIME)
        #vview.set_title(title)
        #vview.show(xvel_prev_newton, yvel_prev_newton, EPS_LOW)
        sprintf(title, "Pressure, time %g", TIME)
        pview.set_title(title)
        pview.show(p_prev)

        # Copy the result of the Newton's iteration into the 
        # previous time level solutions
        xvel_prev_time.copy(xvel_prev_newton)
        yvel_prev_time.copy(yvel_prev_newton)
    else:
        # Assemble and solve
        xvel_sln = Solution()
        yvel_sln = Solution()
        p_sln = Solution()
      
        ls.assemble()
        ls.solve_system(xvel_sln, yvel_sln, p_sln)
        
        # Show the solution at the end of time step
        sprintf(title, "Velocity, time %g", TIME)
        #vview.set_title(title)
        #vview.show(xvel_sln, yvel_sln, EPS_LOW)
        sprintf(title, "Pressure, time %g", TIME)
        pview.set_title(title)
        pview.show(p_sln)

        # This copy destroys xvel_sln and yvel_sln 
        # which are no longer needed
        xvel_prev_time = xvel_sln
        yvel_prev_time = yvel_sln
