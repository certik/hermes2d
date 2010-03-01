#! /usr/bin/env python

#  This is another example that allows you to compare h- and hp-adaptivity from the point of view
#  of both CPU time requirements and discrete problem size, look at the quality of the a-posteriori
#  error estimator used by Hermes (exact error is provided), etc. You can also change
#  the parameter MESH_REGULARITY to see the influence of hanging nodes on the adaptive process.
#  The problem is made harder for adaptive algorithms by increasing the parameter SLOPE.
#
#  PDE: -Laplace u = f
#
#  Known exact solution, see functions fn() and fndd()
#
#  Domain: unit square (0, 1)x(0, 1), see the file square.mesh
#
#  BC:  Dirichlet, given by exact solution

# Import modules
from hermes2d import (Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space,
        WeakForm, Solution, DummySolver, LinSystem, ScalarView, RefSystem,
        H1OrthoHP, set_verbose)
from hermes2d.examples.c22 import set_bc, set_forms

# The following parameters can be changed:
threshold = 0.3
strategy = 0

h_only = False
error_tol = 1
interactive_plotting = False    # should the plot be interactively updated
                                # during the calculation? (slower)
show_mesh = True
show_graph = True

set_verbose(False)

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

sln = Solution()
rsln = Solution()
solver = DummySolver()

view = ScalarView("Solution")
mview = MeshView("Mesh")
graph = []
iter = 0
print "Calculating..."

while 1:
    space.assign_dofs()

    sys = LinSystem(wf, solver)
    sys.set_spaces(space)
    sys.set_pss(pss)
    sys.assemble()
    sys.solve_system(sln)
    dofs = sys.get_matrix().shape[0]
    if interactive_plotting:
        view.show(sln, lib="mayavi", filename="a%02d.png" % iter)
        if show_mesh:
            mview.show(mesh, lib="mpl", method="orders", filename="b%02d.png" % iter)

    rsys = RefSystem(sys)
    rsys.assemble()

    rsys.solve_system(rsln)

    hp = H1OrthoHP(space)
    error_est =  hp.calc_error(sln, rsln)*100
    print "iter=%02d, error_est=%5.2f%%, DOFS=%d" % (iter, error_est, dofs)
    graph.append([dofs, error_est])
    if error_est < error_tol:
        break
    hp.adapt(threshold, strategy, h_only)
    iter += 1

if not interactive_plotting:
    view.show(sln, lib="mayavi")
    if show_mesh:
        mview = MeshView("Mesh")
        mview.show(mesh, lib="mpl", method="orders")
        mview.wait()

if show_graph:
    from numpy import array
    graph = array(graph)
    import pylab
    pylab.clf()
    pylab.plot(graph[:, 0], graph[:, 1], "ko", label="error estimate")
    pylab.plot(graph[:, 0], graph[:, 1], "k-")
    pylab.title("Error Convergence for the Inner Layer Problem")
    pylab.legend()
    pylab.xlabel("Degrees of Freedom")
    pylab.ylabel("Error [%]")
    pylab.yscale("log")
    pylab.grid()
    pylab.savefig("graph.png")
