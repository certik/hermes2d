from hermes2d import Linearizer, Solution

def sln2png(sln, filename):
    """
    Creates a nice png image of the Solution sln.
    """
    plot_sln_mayavi(sln)
    from enthought.mayavi.mlab import savefig
    savefig(filename)

def plot_sln_mpl(sln, method="default", just_mesh=False):
    """
    Plots the Solution() instance sln using Linearizer() and matplotlib.

    method = "default" ... creates a plot using triangles (the triangles are
                not interpolated, so sometimes one can see small defects)
    method = "contour" ... takes the vertices from linearizer and interpolates
                them using contour and contourf (it doesn't take into account
                the triangulation, so one can see the defects from the convex
                hull approximation)

    just_mesh ... only shows the mesh, but not the solution
    """
    lin = Linearizer()
    lin.process_solution(sln)
    v = lin.get_vertices()
    if method=="contour":
        from numpy import min, max, linspace
        from matplotlib.mlab import griddata
        import matplotlib.pyplot as plt
        x = v[:, 0]
        y = v[:, 1]
        z = v[:, 2]
        # define grid.
        xi = linspace(min(x), max(x), 100)
        yi = linspace(min(y), max(y), 100)
        # grid the data.
        zi = griddata(x, y, z, xi, yi)
        # contour the gridded data, plotting dots at the nonuniform data points.
        CS = plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
        CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.jet)
        plt.colorbar()
        plt.title('Solution')
    elif method == "default":
        from numpy import array
        import matplotlib.collections as collections
        import matplotlib.pyplot as plt
        verts = []
        vals = []
        for t in lin.get_triangles():
            triangle = tuple([tuple(v[n][:2]) for n in t])
            val = sum([v[n][2] for n in t])
            vals.append(val/3.)
            verts.append(triangle)
        verts = array(verts)
        vals = array(vals)
        if just_mesh:
            lw = 1
        else:
            lw = 0
        col = collections.PolyCollection(verts, linewidths=lw, antialiaseds=0)
        col.set_array(vals)
        col.set_cmap(plt.cm.jet)
        import pylab
        pylab.ion()
        fig = pylab.gcf()#plt.figure()
        pylab.ioff()
        ax = fig.gca()
        ax.add_collection(col)
        ax.set_xlim(verts[:, :, 0].min(), verts[:, :, 0].max())
        ax.set_ylim(verts[:, :, 1].min(), verts[:, :, 1].max())
        ax.set_aspect("equal")
        #plt.colorbar()
        plt.title('Solution')
    else:
        raise ValueError("Unknown method (%s)" % method)

def plot_sln_mayavi(sln, notebook=False):
    """
    Plots the Solution() instance sln using Linearizer() and matplotlib.

    Currently only a very simple version is implemented, that takes the
    vertices from linearizer and interpolates them. More sophisticated version
    should take the triangles.
    """
    lin = Linearizer()
    lin.process_solution(sln)
    vert = lin.get_vertices()
    triangles = lin.get_triangles()
    from numpy import zeros
    from enthought.mayavi import mlab
    x = vert[:, 0]
    y = vert[:, 1]
    z = zeros(len(y))
    t = vert[:, 2]
    if notebook:
        # the off screen rendering properly works only with VTK-5.2 or above:
        mlab.options.offscreen = True
    s = mlab.triangular_mesh(x, y, z, triangles, scalars=t)
    mlab.view(0, 0)

    # Below is a code that does exactly what the "View along the +Z axis"
    # button does:
    #scene = mlab.get_engine().current_scene.scene
    #scene.camera.focal_point = [0, 0, 0]
    #scene.camera.position = [0, 0, 1]
    #scene.camera.view_up = [0, 1, 0]
    #scene.renderer.reset_camera()
    #scene.render()
    # the above looks ok, but there is still quite a large margin, so we prefer
    # to just call .view(0, 0), which seems to be working fine.
    return s

def plot_hermes_mesh_mpl(mesh):
    # XXX: extract the mesh from hermes here:
    nodes = ((0.6, 0), (0.7, 0), (0.59396969619669993, 0.084852813742385721),
            (0.7, 0.1), (0.51428571428571435, 0.30904725218262763), (0.6,
                0.36055512754639885), (0.7, 0.36055512754639885), (0, 0.6), (0,
                    0.7), (0.084852813742385721, 0.59396969619669993), (0.1,
                        0.7), (0.30904725218262763, 0.51428571428571435),
                    (0.36055512754639885, 0.6), (0.36055512754639885, 0.7),
                    (0.6, 0.6), (0.7, 0.6), (0.7, 0.7), (0.6, 0.7), (0.7,
                        -0.1), (0.62928932188134523, -0.070710678118654766),
                    (-0.1, 0.7), (-0.070710678118654766, 0.62928932188134523))

    elements = ((0, 1, 3, 2), (2, 3, 5, 4), (6, 5, 3), (8, 7, 9, 10), (10, 9,
        11, 12), (13, 10, 12), (4, 5, 12, 11), (5, 6, 15, 14), (13, 12, 14,
            17), (14, 15, 16, 17), (0, 19, 1), (19, 18, 1), (21, 7, 8), (20,
                21, 8))

    return plot_mesh_mpl(nodes, elements)

def plot_mesh_mpl(nodes, elements, polynomial_orders=None, colors=None):
    import numpy as np
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    import matplotlib.pyplot as pyplot

    if polynomial_orders is None:
        polynomial_orders = [1] * len(elements)

    if colors is None:
        colors = {0: '#000000', 1: '#000684', 2: '#3250fc',
            3: '#36c4ee', 4: '#04eabc', 5: '#62ff2a', 6: '#fdff07',
            7: '#ffa044', 8: '#ff1111', 9: '#b02c2c', 10: '#820f97'}

    # check that if orders and elements match (if orders are passed in)
    if polynomial_orders is not None:
        assert len(elements) == len(polynomial_orders)

    path_polynomial_orders = {}
    pathpatch_polynomial_orders = {}
    vertices_polynomial_orders = {}
    codes_polynomial_orders = {}

    for key, value in colors.items():
        if not key in polynomial_orders:
            continue
        path_polynomial_orders[key] = 'path will be added later'
        pathpatch_polynomial_orders[key] = 'patchPath will be added later'
        vertices_polynomial_orders[key] = []
        codes_polynomial_orders[key] = []

    # join nodes with lines:
    for i, e in enumerate(elements):
        x_avg = 0
        y_avg = 0
        for k, node_index in enumerate(e):
            vertices_polynomial_orders[polynomial_orders[i]].append(nodes[node_index])
            if k == 0:
                codes_polynomial_orders[polynomial_orders[i]].append(Path.MOVETO)
            else:
                codes_polynomial_orders[polynomial_orders[i]].append(Path.LINETO)

        vertices_polynomial_orders[polynomial_orders[i]].append((0,0))
        codes_polynomial_orders[polynomial_orders[i]].append(Path.CLOSEPOLY)


    for key, color in colors.items():
        if not key in polynomial_orders:
            continue
        vertices_polynomial_orders[key] = np.array(vertices_polynomial_orders[key], float)
        path_polynomial_orders[key] = Path(vertices_polynomial_orders[key], codes_polynomial_orders[key])
        pathpatch_polynomial_orders[key] = PathPatch(path_polynomial_orders[key], facecolor=color, edgecolor='green')


    fig = pyplot.figure()
    sp = fig.add_subplot(111)

    for key, patch in pathpatch_polynomial_orders.items():
        sp.add_patch(patch)

    for i, e in enumerate(elements):
        x_avg = 0
        y_avg = 0
        for k, node_index in enumerate(e):
            x1, y1 = nodes[e[k-1]]
            x2, y2 = nodes[e[k]]

            x_avg += x2
            y_avg += y2

        x_avg /= len(e)
        y_avg /= len(e)
        sp.text(x_avg, y_avg, str(polynomial_orders[i]))

    #for key,vertices in vertices_polynomial_orders.items():
    #    sp.dataLim.update_from_data_xy(vertices)

    sp.set_title('Mesh using path & Patches')

    sp.autoscale_view()
    return sp.figure


class ScalarView(object):

    def __init__(self, x=0, y=0, w=50, h=50, name="Solution"):
        self._name = name
        self._lib = None
        self._notebook = False

    def show_scale(self, *args):
        pass

    def show_mesh(self, *args):
        pass

    def wait(self):
        if self._lib == "mpl" and self._notebook == False:
            import pylab
            pylab.show()

    def show(self, sln, show=True, lib="mpl", notebook=False, filename="a.png",
            **options):
        """
        Shows the solution.

        show ... should it actually plot the window? Set to False in tests.
        lib .... which library to use for the plotting? either "mpl" or "mayavi"
        notebook ... are we running inside Sage notebook? If True, just save
                the image to a.png
        filename ... the name of the filename if we are saving the image (e.g.
                notebook == False)
        """
        self._lib = lib
        self._notebook = notebook
        if lib == "mpl":
            plot_sln_mpl(sln, **options)
            import pylab
            if show:
                if notebook:
                    pylab.savefig(filename)
                else:
                    pylab.ion()
                    pylab.draw()
                    pylab.ioff()
        elif lib == "mayavi":
            plot_sln_mayavi(sln, notebook=notebook)
            from enthought.mayavi import mlab
            if show:
                if notebook:
                    mlab.savefig(filename)
                else:
                    mlab.show()
        else:
            raise NotImplementedError("Unknown library '%s'" % lib)

class MeshView(object):

    def __init__(self, name="Solution", x=0, y=0, w=500, h=500):
        self._name = name
        self._x = x
        self._y = y
        self._w = w
        self._h = h

    def wait(self):
        pass

    def show(self, mesh, show=True, lib="glut", notebook=False,
            filename="a.png", **options):
        if lib == "glut":
            from _hermes2d import MeshView
            m = MeshView(self._name, self._x, self._y, self._w, self._h)
            m.show(mesh)
            m.wait()
        elif lib == "mpl":
            p = plot_hermes_mesh_mpl(mesh)
            if show:
                if notebook:
                    p.savefig(filename)
                else:
                    p.show()
            return p
        else:
            raise NotImplementedError("Unknown library '%s'" % lib)
