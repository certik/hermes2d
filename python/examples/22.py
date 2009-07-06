from hermes2d import (Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space,
        WeakForm, Solution, DummySolver, LinSystem, ScalarView, RefSystem)

from hermes2d.examples.c22 import set_bc, set_forms

# Create a mesh:
mesh = Mesh()
mesh.create([
        [0, 0],
        [1, 0],
        [1, 1],
        [0, 1],
    ], [
        [2, 3, 0, 1, 0],
    ], [
        [0, 1, 1],
        [1, 2, 1],
        [2, 3, 1],
        [3, 0, 1],
    ], [])

mesh.refine_all_elements()

shapeset = H1Shapeset()
pss = PrecalcShapeset(shapeset)

space = H1Space(mesh, shapeset)
set_bc(space)
space.set_uniform_order(1)

wf = WeakForm(1)
set_forms(wf)

space.assign_dofs()

sln = Solution()
solver = DummySolver()
sys = LinSystem(wf, solver)
sys.set_spaces(space)
sys.set_pss(pss)
sys.assemble()
sys.solve_system(sln)

rsln = Solution()
rsys = RefSystem(sys)
rsys.assemble()

#A = rsys.get_matrix()
#rhs = rsys.get_rhs()
#from scipy.sparse.linalg import cg
#x, res = cg(A, rhs)
#rsln.set_fe_solution(rsys.spaces[0], rsys.pss[0], x)

rsys.solve_system(rsln)

#view = ScalarView("Solution")
#view.show(sln)
