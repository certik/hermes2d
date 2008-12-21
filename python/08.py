#! /usr/bin/env python

from hermes2d import finalize, Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
        DiscreteProblem, Solution, ScalarView, VectorView

from c08 import set_bc, set_forms

mesh = Mesh()
mesh.load("cylinder4.mesh")
#mesh.refine_element(0)
#mesh.refine_all_elements()
mesh.refine_towards_boundary(5, 3)
shapeset = H1Shapeset()
pss = PrecalcShapeset(shapeset)

# create an H1 space
xvel = H1Space(mesh, shapeset)
yvel = H1Space(mesh, shapeset)
press = H1Space(mesh, shapeset)
xvel.set_uniform_order(2)
yvel.set_uniform_order(2)
press.set_uniform_order(1)

set_bc(xvel, yvel, press)

ndofs = 0
ndofs += xvel.assign_dofs(ndofs)
ndofs += yvel.assign_dofs(ndofs)
ndofs += press.assign_dofs(ndofs)

xprev = Solution()
yprev = Solution()

xprev.set_zero(mesh)
yprev.set_zero(mesh)

# initialize the discrete problem
dp = DiscreteProblem()
dp.set_num_equations(3)
dp.set_spaces(xvel, yvel, press)
dp.set_pss(pss)
dp.set_external_fns(xprev, yprev)
set_forms(dp, xprev, yprev)

# visualize the solution
vview = VectorView("velocity [m/s]", 0, 0, 1200, 350)
pview = ScalarView("pressure [Pa]", 0, 500, 1200, 350)
vview.set_min_max_range(0, 1.9)
vview.show_scale(False)
pview.show_scale(False)
pview.show_mesh(False)

# assemble the stiffness matrix and solve the system
dp.create_matrix();

EPS_LOW = 0.0014

for i in range(1000):
    print "*** Iteration %d ***" % i
    psln = Solution()
    dp.assemble_matrix_and_rhs()
    dp.solve_system(xprev, yprev, psln)
    vview.show(xprev, yprev, 2*EPS_LOW)
    pview.show(psln)

finalize()
