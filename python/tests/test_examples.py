def test_example_01():
    from hermes2d import finalize, Mesh, MeshView
    mesh = Mesh()
    mesh.load("domain.mesh")
    mesh.refine_all_elements()
    mesh.refine_all_elements()
    mesh.refine_all_elements()

def test_example_02():
    from hermes2d import finalize, Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
            DiscreteProblem, BaseView, set_verbose

    from c02 import set_forms

    set_verbose(False)

    mesh = Mesh()
    mesh.load("domain.mesh")
    mesh.refine_element(0)
    shapeset = H1Shapeset()
    pss = PrecalcShapeset(shapeset)

    # create an H1 space
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(5)
    space.assign_dofs();

    # initialize the discrete problem
    dp = DiscreteProblem()
    dp.set_num_equations(1)
    dp.set_spaces(space)
    dp.set_pss(pss)
    set_forms(dp)

def test_example_03():
    from hermes2d import finalize, Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
            DiscreteProblem, Solution, ScalarView, set_verbose

    from c02 import set_forms

    set_verbose(False)

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

def test_example_04():
    from hermes2d import finalize, Mesh, H1Shapeset, PrecalcShapeset, \
            H1Space, DiscreteProblem, Solution, ScalarView, set_verbose

    from c04 import set_bc, set_forms

    set_verbose(False)

    mesh = Mesh()
    mesh.load("domain.mesh")
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
    dp = DiscreteProblem()
    dp.set_num_equations(1)
    dp.set_spaces(space)
    dp.set_pss(pss)
    set_forms(dp)

    sln = Solution()
    dp.create_matrix()
    dp.assemble_matrix_and_rhs()
    dp.solve_system(sln)

def test_example_05():
    from hermes2d import finalize, Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
            DiscreteProblem, Solution, ScalarView, set_verbose

    from c05 import set_bc, set_forms

    set_verbose(False)

    mesh = Mesh()
    mesh.load("domain.mesh")
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
    dp = DiscreteProblem()
    dp.set_num_equations(1)
    dp.set_spaces(space)
    dp.set_pss(pss)
    set_forms(dp)

    sln = Solution()
    dp.create_matrix()
    dp.assemble_matrix_and_rhs()
    dp.solve_system(sln)

def test_example_06():
    from hermes2d import finalize, Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
            DiscreteProblem, Solution, ScalarView, set_verbose

    from c06 import set_bc, set_forms

    set_verbose(False)

    mesh = Mesh()
    mesh.load("domain.mesh")
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
    dp = DiscreteProblem()
    dp.set_num_equations(1)
    dp.set_spaces(space)
    dp.set_pss(pss)
    set_forms(dp)

    sln = Solution()
    dp.create_matrix()
    dp.assemble_matrix_and_rhs()
    dp.solve_system(sln)

def test_example_07():
    from hermes2d import finalize, Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
            DiscreteProblem, Solution, ScalarView, VonMisesFilter, \
            set_verbose, set_warn_integration

    from c07 import set_bc, set_forms

    set_verbose(False)

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
    old_flag = set_warn_integration(False)
    dp.assemble_matrix_and_rhs()
    set_warn_integration(old_flag)
    dp.solve_system(xsln, ysln)

    E = float(200e9)
    nu = 0.3
    stress = VonMisesFilter(xsln, ysln, E / (2*(1 + nu)),
            (E * nu) / ((1 + nu) * (1 - 2*nu)))


def test_example_08():
    from hermes2d import finalize, Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
            DiscreteProblem, Solution, ScalarView, VectorView, set_verbose

    from c08 import set_bc, set_forms

    set_verbose(False)

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

    # assemble the stiffness matrix and solve the system
    dp.create_matrix();

    EPS_LOW = 0.0014

    for i in range(3):
        psln = Solution()
        dp.assemble_matrix_and_rhs()
        dp.solve_system(xprev, yprev, psln)
