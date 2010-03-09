from hermes2d import Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
        LinSystem, WeakForm, Solution, ScalarView, VonMisesFilter, \
        set_verbose, set_warn_integration, DummySolver, EPS_HIGH, FN_DX, \
        FN_DY, H1OrthoHP, RefSystem
from hermes2d.forms import set_forms
from hermes2d.examples import get_example_mesh, get_sample_mesh, \
        get_cylinder_mesh, get_07_mesh, get_cathedral_mesh, get_bracket_mesh

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
    from hermes2d.examples.c03 import set_bc
    set_bc(space)
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

    P_INIT = 2             # Initial polynomial degree of all mesh elements.

    mesh = Mesh()
    mesh.load(get_07_mesh())

    # Initialize the shapeset and the cache
    shapeset = H1Shapeset()
    pss = PrecalcShapeset(shapeset)

    #create finite element space
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(P_INIT)
    set_bc(space)

    # Enumerate basis functions
    space.assign_dofs()

    # weak formulation
    wf = WeakForm(1)
    set_forms(wf)

    #matrix solver
    solver = DummySolver()

    #Solve the problem
    sln = Solution()
    ls = LinSystem(wf, solver)
    ls.set_spaces(space)
    ls.set_pss(pss)
    ls.assemble()
    ls.solve_system(sln)


def test_example_08():
    from hermes2d.examples.c08 import set_bc, set_forms

    set_verbose(False)

    mesh = Mesh()
    mesh.load(get_sample_mesh())

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
    sys.assemble()
    sys.solve_system(xsln, ysln)

    E = float(200e9)
    nu = 0.3
    stress = VonMisesFilter(xsln, ysln, E / (2*(1 + nu)), (E * nu) / ((1 + nu) * (1 - 2*nu)))

def test_example_09():
    from hermes2d.examples.c09 import set_bc, temp_ext, set_forms

    # The following parameters can be played with:
    P_INIT = 1            # polynomial degree of elements
    INIT_REF_NUM = 4      # number of initial uniform refinements
    TAU = 300.0           # time step in seconds

    # Problem constants
    T_INIT = 10           # temperature of the ground (also initial temperature)
    FINAL_TIME = 86400    #length of time interval (24 hours) in seconds

    # Global variable
    TIME = 0;

    # Load the mesh
    mesh = Mesh()
    mesh.load(get_cathedral_mesh())

    #for i in range(INIT_REF_NUM):
    #    mesh.refine_all_elements()
    #mesh.refine_towards_boundary(2, 5)

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

    # Assemble and solve
    ls.assemble()
    rhsonly = True
    ls.solve_system(tsln)

def test_example_10():
    from hermes2d.examples.c10 import set_bc, set_forms
    from hermes2d.examples import get_motor_mesh

    # The following parameters can be changed:

    P_INIT = 1              # Initial polynomial degree of all mesh elements.
    THRESHOLD = 0.2         # This is a quantitative parameter of the adapt(...) function and
                            # it has different meanings for various adaptive strategies (see below).

    STRATEGY = 1            # Adaptive strategy:
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

    ERR_STOP = 0.1          # Stopping criterion for adaptivity (rel. error tolerance between the
                            # fine mesh and coarse mesh solution in percent).

    NDOF_STOP = 40000       # Adaptivity process stops when the number of degrees of freedom grows
                            # over this limit. This is to prevent h-adaptivity to go on forever.

    # Load the mesh
    mesh = Mesh()
    mesh.load(get_motor_mesh())

    # Initialize the shapeset and the cache
    shapeset = H1Shapeset()
    pss = PrecalcShapeset(shapeset)

    # Create finite element space
    space = H1Space(mesh, shapeset)
    set_bc(space)
    space.set_uniform_order(P_INIT)

    # Enumerate basis functions
    space.assign_dofs()

    # Initialize the discrete problem
    wf = WeakForm(1)
    set_forms(wf)

    # Matrix solver
    solver = DummySolver()

    # Adaptivity loop
    it = 1
    ndofs = 0

    done = False
    cpu = 0.0

    sln_coarse = Solution()
    sln_fine = Solution()

    # Solve the coarse mesh problem
    ls = LinSystem(wf, solver)
    ls.set_spaces(space)
    ls.set_pss(pss)
    ls.assemble()
    ls.solve_system(sln_coarse)

    # Solve the fine mesh problem
    rs = RefSystem(ls)
    rs.assemble()
    rs.solve_system(sln_fine)

    # Calculate element errors and total error estimate
    hp = H1OrthoHP(space);
    err_est = hp.calc_error(sln_coarse, sln_fine) * 100

def test_example_11():
    from hermes2d.examples.c11 import set_bc, set_wf_forms, set_hp_forms

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
    ADAPT_TYPE = 0           # Type of automatic adaptivity:
                                # ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                # ADAPT_TYPE = 1 ... adaptive h-FEM,
                                # ADAPT_TYPE = 2 ... adaptive p-FEM.
    ISO_ONLY = False         # Isotropic refinement flag (concerns quadrilateral elements only).
                                # ISO_ONLY = false ... anisotropic refinement of quad elements
                                # is allowed (default),
                                # ISO_ONLY = true ... only isotropic refinements of quad elements
                                # are allowed.
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

    # initial mesh refinements
    xmesh.refine_element(1)
    xmesh.refine_element(4)

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

    # Calculating the number of degrees of freedom
    ndofs = xdisp.assign_dofs()
    ndofs += ydisp.assign_dofs(ndofs)
    
    # Solve the coarse mesh problem
    ls = LinSystem(wf, solver)
    ls.set_spaces(xdisp, ydisp)
    ls.set_pss(xpss, ypss)
    ls.assemble()
    ls.solve_system(x_sln_coarse, y_sln_coarse)

    # View the solution -- this can be slow; for illustration only
    stress_coarse = VonMisesFilter(x_sln_coarse, y_sln_coarse, mu, lamda)

    # Solve the fine mesh problem
    rs = RefSystem(ls)
    rs.assemble()
    rs.solve_system(x_sln_fine, y_sln_fine)

    # Calculate element errors and total error estimate
    hp = H1OrthoHP(xdisp, ydisp)
    set_hp_forms(hp)
    err_est = hp.calc_error_2(x_sln_coarse, y_sln_coarse, x_sln_fine, y_sln_fine) * 100

    # Show the fine solution - this is the final result
    stress_fine = VonMisesFilter(x_sln_fine, y_sln_fine, mu, lamda)

def test_example_12():
    from hermes2d.examples.c12 import set_bc, set_forms
    from hermes2d.examples import get_example_mesh

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
    mesh.load(get_example_mesh())
    #mesh.load("hermes2d/examples/12.mesh")

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

    # Matrix solver
    solver = DummySolver()

    # Adaptivity loop
    it = 0
    ndofs = 0
    done = False
    sln_coarse = Solution()
    sln_fine = Solution()

    # Solve the coarse mesh problem
    ls = LinSystem(wf, solver)
    ls.set_spaces(space)
    ls.set_pss(pss)
    
    ls.assemble()
    ls.solve_system(sln_coarse)

    # Solve the fine mesh problem
    rs = RefSystem(ls)
    rs.assemble()
    rs.solve_system(sln_fine)
    
    # Calculate element errors and total error estimate
    hp = H1OrthoHP(space);
    err_est = hp.calc_error(sln_coarse, sln_fine) * 100
