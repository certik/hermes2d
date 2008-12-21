#! /usr/bin/env python

from hermes2d import finalize, Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
        DiscreteProblem, Solution, ScalarView, VonMisesFilter

from c07 import set_bc, set_forms

mesh = Mesh()
mesh.load("sample.mesh")
#mesh.refine_element(0)
#mesh.refine_all_elements()
#mesh.refine_towards_boundary(5, 3)
shapeset = H1Shapeset()
pss = PrecalcShapeset(shapeset)

# create an H1 space
xdisp = H1Space(mesh, shapeset)
ydisp = H1Space(mesh, shapeset)
xdisp.set_uniform_order(8)
ydisp.set_uniform_order(8)

set_bc(xdisp, ydisp)

ndofs = xdisp.assign_dofs(0)
ndofs += ydisp.assign_dofs(ndofs)

# initialize the discrete problem
dp = DiscreteProblem()
dp.set_num_equations(2)
dp.set_spaces(xdisp, ydisp)
dp.set_pss(pss)
set_forms(dp)

xsln = Solution()
ysln = Solution()
dp.create_matrix()
dp.assemble_matrix_and_rhs()
dp.solve_system(xsln, ysln)

view = ScalarView("Von Mises stress [Pa]", 50, 50, 1200, 600)
E = float(200e9)
nu = 0.3
stress = VonMisesFilter(xsln, ysln, E / (2*(1 + nu)),
        (E * nu) / ((1 + nu) * (1 - 2*nu)))
view.show(stress)

finalize()
