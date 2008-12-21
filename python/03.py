#! /usr/bin/env python

from hermes2d import finalize, Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
        DiscreteProblem, Solution, ScalarView

from c02 import set_forms

mesh = Mesh()
mesh.load("domain.mesh")
mesh.refine_element(0)
shapeset = H1Shapeset()
pss = PrecalcShapeset(shapeset)

# create an H1 space
space = H1Space(mesh, shapeset)
space.set_uniform_order(5)
space.assign_dofs()

# initialize the discrete problem
dp = DiscreteProblem()
dp.set_num_equations(1)
dp.set_spaces(space)
dp.set_pss(pss)
set_forms(dp)

# assemble the stiffness matrix and solve the system
sln = Solution()
dp.create_matrix()
dp.assemble_matrix_and_rhs()
dp.solve_system(sln)

# visualize the solution
view = ScalarView("Solution")
view.show(sln)

finalize()
