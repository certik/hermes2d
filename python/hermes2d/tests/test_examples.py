from hermes2d import Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
        LinSystem, WeakForm, Solution, ScalarView, VonMisesFilter, \
        set_verbose, set_warn_integration, DummySolver
from hermes2d.forms import set_forms
from hermes2d.examples import get_example_mesh, get_sample_mesh, \
        get_cylinder_mesh

domain_mesh = get_example_mesh()
sample_mesh = get_sample_mesh()
cylinder_mesh = get_cylinder_mesh()

def test_example_01():
    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_all_elements()
    mesh.refine_all_elements()
    mesh.refine_all_elements()

def test_example_02():
    set_verbose(False)

    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_element(0)
    shapeset = H1Shapeset()
    pss = PrecalcShapeset(shapeset)

    # create an H1 space
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(5)
    space.assign_dofs();

    # initialize the discrete problem
    wf = WeakForm(1)
    set_forms(wf)

    solver = DummySolver()
    sys = LinSystem(wf, solver)
    sys.set_spaces(space)
    sys.set_pss(pss)

def test_example_03():
    set_verbose(False)

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
    sln = Solution()
    sys.assemble()
    sys.solve_system(sln)
    assert abs(sln.l2_norm() - 0.25493) < 1e-4
    assert abs(sln.h1_norm() - 0.89534) < 1e-4

def test_example_04():
    from hermes2d.examples.c04 import set_bc

    set_verbose(False)

    mesh = Mesh()
    mesh.load(domain_mesh)
    #mesh.refine_element(0)
    #mesh.refine_all_elements()
    mesh.refine_towards_boundary(5, 3)
    shapeset = H1Shapeset()
    pss = PrecalcShapeset(shapeset)

    # create an H1 space
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(5)

    set_bc(space)

    space.assign_dofs()

    xprev = Solution()
    yprev = Solution()

    # initialize the discrete problem
    wf = WeakForm()
    set_forms(wf, -4)

    solver = DummySolver()
    sys = LinSystem(wf, solver)
    sys.set_spaces(space)
    sys.set_pss(pss)

    # assemble the stiffness matrix and solve the system
    sys.assemble()
    sln = Solution()
    sys.solve_system(sln)
    assert abs(sln.l2_norm() - 1.22729) < 1e-4
    assert abs(sln.h1_norm() - 2.90006) < 1e-4

def test_example_05():
    from hermes2d.examples.c05 import set_bc
    from hermes2d.examples.c05 import set_forms as set_forms_surf

    set_verbose(False)

    mesh = Mesh()
    mesh.load(domain_mesh)
    mesh.refine_towards_vertex(3, 12)
    shapeset = H1Shapeset()
    pss = PrecalcShapeset(shapeset)

    # create an H1 space
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(4)

    set_bc(space)

    space.assign_dofs()

    xprev = Solution()
    yprev = Solution()

    # initialize the discrete problem
    wf = WeakForm(1)
    set_forms(wf, -1)
    set_forms_surf(wf)

    sln = Solution()
    solver = DummySolver()
    sys = LinSystem(wf, solver)
    sys.set_spaces(space)
    sys.set_pss(pss)
    sys.assemble()
    sys.solve_system(sln)
    assert abs(sln.l2_norm() - 0.535833) < 1e-4
    assert abs(sln.h1_norm() - 1.332908) < 1e-4

def test_example_06():
    from hermes2d.examples.c06 import set_bc, set_forms

    set_verbose(False)

    mesh = Mesh()
    mesh.load(domain_mesh)
    #mesh.refine_element(0)
    #mesh.refine_all_elements()
    mesh.refine_towards_boundary(5, 3)
    shapeset = H1Shapeset()
    pss = PrecalcShapeset(shapeset)

    # create an H1 space
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(5)

    set_bc(space)

    space.assign_dofs()

    xprev = Solution()
    yprev = Solution()

    # initialize the discrete problem
    wf = WeakForm(1)
    set_forms(wf)

    solver = DummySolver()
    sys = LinSystem(wf, solver)
    sys.set_spaces(space)
    sys.set_pss(pss)

    sln = Solution()
    sys.assemble()
    sys.solve_system(sln)
    assert abs(sln.l2_norm() - 121.78788) < 1e-4
    assert abs(sln.h1_norm() - 126.96528) < 1e-4

def test_example_07():
    from hermes2d.examples.c07 import set_bc, set_forms

    set_verbose(False)

    mesh = Mesh()
    mesh.load(sample_mesh)
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
    wf = WeakForm(2)
    set_forms(wf)

    solver = DummySolver()
    sys = LinSystem(wf, solver)
    sys.set_spaces(xdisp, ydisp)
    sys.set_pss(pss)

    xsln = Solution()
    ysln = Solution()
    old_flag = set_warn_integration(False)
    sys.assemble()
    set_warn_integration(old_flag)
    sys.solve_system(xsln, ysln)

    E = float(200e9)
    nu = 0.3
    stress = VonMisesFilter(xsln, ysln, E / (2*(1 + nu)),
            (E * nu) / ((1 + nu) * (1 - 2*nu)))


def test_example_08():
    from hermes2d.examples.c08 import set_bc, set_forms

    set_verbose(False)

    mesh = Mesh()
    mesh.load(cylinder_mesh)
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
    wf = WeakForm(3)
    set_forms(wf, xprev, yprev)

    solver = DummySolver()
    sys = LinSystem(wf, solver)
    sys.set_spaces(xvel, yvel, press)
    sys.set_pss(pss)
    #dp.set_external_fns(xprev, yprev)

    # visualize the solution

    EPS_LOW = 0.0014

    for i in range(3):
        psln = Solution()
        sys.assemble()
        sys.solve_system(xprev, yprev, psln)
