from hermes2d import Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
        WeakForm, Solution, ScalarView, LinSystem, DummySolver, raises

from hermes2d.forms import set_forms

from hermes2d.examples import get_example_mesh

domain_mesh = get_example_mesh()

def test_ScalarView_mpl_default():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_element(0)
    shapeset = H1Shapeset()
    pss = PrecalcShapeset(shapeset)

    # create an H1 space
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(5)
    space.assign_dofs()

    # initialize the discrete problem
    wf = WeakForm(1)
    set_forms(wf)

    solver = DummySolver()
    sys = LinSystem(wf, solver)
    sys.set_spaces(space)
    sys.set_pss(pss)

    # assemble the stiffness matrix and solve the system
    sys.assemble()
    A = sys.get_matrix()
    b = sys.get_rhs()
    from scipy.sparse.linalg import cg
    x, res = cg(A, b)
    sln = Solution()
    sln.set_fe_solution(space, pss, x)

    view = ScalarView("Solution")
    view.show(sln, show=False)

def test_ScalarView_mpl_default():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_element(0)
    shapeset = H1Shapeset()
    pss = PrecalcShapeset(shapeset)

    # create an H1 space
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(5)
    space.assign_dofs()

    # initialize the discrete problem
    wf = WeakForm(1)
    set_forms(wf)

    solver = DummySolver()
    sys = LinSystem(wf, solver)
    sys.set_spaces(space)
    sys.set_pss(pss)

    # assemble the stiffness matrix and solve the system
    sys.assemble()
    A = sys.get_matrix()
    b = sys.get_rhs()
    from scipy.sparse.linalg import cg
    x, res = cg(A, b)
    sln = Solution()
    sln.set_fe_solution(space, pss, x)

    view = ScalarView("Solution")
    view.show(sln, show=False, method="contour")

def test_ScalarView_mpl_unknown():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_element(0)
    shapeset = H1Shapeset()
    pss = PrecalcShapeset(shapeset)

    # create an H1 space
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(5)
    space.assign_dofs()

    # initialize the discrete problem
    wf = WeakForm(1)
    set_forms(wf)

    solver = DummySolver()
    sys = LinSystem(wf, solver)
    sys.set_spaces(space)
    sys.set_pss(pss)

    # assemble the stiffness matrix and solve the system
    sys.assemble()
    A = sys.get_matrix()
    b = sys.get_rhs()
    from scipy.sparse.linalg import cg
    x, res = cg(A, b)
    sln = Solution()
    sln.set_fe_solution(space, pss, x)

    view = ScalarView("Solution")
    assert raises(ValueError, 'view.show(sln, show=False, method="something_unknown_123")')
