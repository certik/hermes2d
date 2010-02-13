from hermes2d import Linearizer, Solution

def sln2png(sln, filename):
    """
    Creates a nice png image of the Solution sln.
    """
    plot_sln_mayavi(sln)
    from enthought.mayavi.mlab import savefig
    savefig(filename)

def plot_sln_mpl(sln, method="default", just_mesh=False, axes=None):
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
        #import matplotlib.pyplot as plt
        if axes is None:
            from pylab import gca
            axes = gca()
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
        #col.set_cmap(plt.cm.jet)
        ax = axes
        ax.add_collection(col)
        ax.set_xlim(verts[:, :, 0].min(), verts[:, :, 0].max())
        ax.set_ylim(verts[:, :, 1].min(), verts[:, :, 1].max())
        ax.set_aspect("equal")
        #plt.colorbar()
        #plt.title('Solution')
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

def plot_hermes_mesh_mpl(mesh, space=None, method="simple"):
    if method == "simple":
        return plot_mesh_mpl_simple(mesh.nodes_dict, mesh.elements)
    elif method == "orders":
        if space != None:
            return plot_mesh_mpl_orders(mesh.nodes_dict, mesh.elements, mesh.curves, polynomial_orders=mesh.get_elements_order(space))
        else:
            return plot_mesh_mpl_orders(mesh.nodes_dict, mesh.elements, mesh.curves)
    else:
        raise ValueError("Unknown method")

def plot_mesh_mpl_orders(nodes, elements, curves=None, polynomial_orders=None, colors=None):
    """
    This plots the mesh together with polynomial orders.

    >>> a = 5
    >>> b = 3
    >>> f = plot_mesh_mpl_orders([(0,-a), (a,-a), (-a,0), (0,0), (a,0), (-a,a), (0,a), (a*b,a*b) ], \
    [ (0, 1, 4, 3), (3, 4, 7), (3, 7, 6), (2, 3, 6, 5) ], None, [1,2,3,4])
    >>> f
    <matplotlib.figure.Figure object at 0x...>

    """
    import math

    import numpy as np
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D

    import matplotlib.pyplot as pyplot

    if colors is None:
        colors = {1: '#000684', 2: '#3250fc',
            3: '#36c4ee', 4: '#04eabc', 5: '#62ff2a', 6: '#fdff07',
            7: '#ffa044', 8: '#ff1111', 9: '#b02c2c', 10: '#820f97'}

    # Check that if orders and elements match (if orders are passed in).
    if polynomial_orders is not None:
        assert len(elements) == len(polynomial_orders)
    else:
        polynomial_orders = [1] * len(elements)

    fig = pyplot.figure()
    sp = fig.add_subplot(111)

    path_polynomial_orders = {}
    pathpatch_polynomial_orders = {}
    vertices_polynomial_orders = {}
    codes_polynomial_orders = {}

    for key,value in colors.items():
        path_polynomial_orders[key] = None
        pathpatch_polynomial_orders[key] = None
        vertices_polynomial_orders[key] = []
        codes_polynomial_orders[key] = []

    curvilinear_elements_index = []
    if curves != None:
        curvilinear_elements_index = [k for k, v in curves.items()]

    # Plot non-curvilinear elements
    for i, e in enumerate(elements):
        # Plot only if this element do not contain curves
        if i not in curvilinear_elements_index:
            for k, node_index in enumerate(e):
                j = polynomial_orders[i]
                vertices_polynomial_orders[j].append(nodes[node_index])
                if k == 0:
                    codes_polynomial_orders[j].append(Path.MOVETO)
                else:
                    codes_polynomial_orders[j].append(Path.LINETO)
            vertices_polynomial_orders[j].append((0, 0))
            codes_polynomial_orders[j].append(Path.CLOSEPOLY)

    # Plot curvilinear elements
    for cei in curvilinear_elements_index:
        for i, e in enumerate(elements):
            # Plot only if this element contains curves
            if i == cei:
                j = polynomial_orders[i]
                vertices_polynomial_orders[j].append(nodes[e[0]])
                codes_polynomial_orders[j].append(Path.MOVETO)
                z = range(len(e)-1)
                z.append(-1)

                for k in z:
                    for cl in curves[cei]:
                        if cl[0] == nodes[e[k]] and cl[2] == nodes[e[k+1]]:
                            vertices_polynomial_orders[j].append(cl[1])
                            codes_polynomial_orders[j].append(Path.CURVE3)
                            vertices_polynomial_orders[j].append(cl[2])
                            codes_polynomial_orders[j].append(Path.CURVE3)
                        else:
                            vertices_polynomial_orders[j].append(nodes[e[k+1]])
                            codes_polynomial_orders[j].append(Path.LINETO)
                vertices_polynomial_orders[j].append((0, 0))
                codes_polynomial_orders[j].append(Path.CLOSEPOLY)

    for key, color in colors.items():
        if len(vertices_polynomial_orders[key]) != 0:
            vertices_polynomial_orders[key] = np.array(vertices_polynomial_orders[key], float)
            path_polynomial_orders[key] = Path(vertices_polynomial_orders[key], codes_polynomial_orders[key])
            pathpatch_polynomial_orders[key] = PathPatch(path_polynomial_orders[key], facecolor=color, edgecolor='#000000')

    for key,patch in pathpatch_polynomial_orders.items():
        if patch != None:
            sp.add_patch(patch)

    # Create legend
    def split_nodes():
        x = []
        y = []

        if isinstance(nodes, dict):
            _nodes = nodes.items()
        else:
            _nodes = enumerate(nodes)
        for k, pnt in _nodes:
            x.append(pnt[0])
            y.append(pnt[1])

        return (x, y)

    def get_max(what='x'):
        x, y = split_nodes()

        if what == 'x':
            return max(x)
        else:
            return max(y)

    def get_min(what='x'):
        x, y = split_nodes()

        if what == 'x':
            return min(x)
        else:
            return min(y)

    maxX = get_max('x')
    maxY = get_max('y')

    minX = get_min('x')
    minY = get_min('y')

    dy = (maxY - minY) / 20
    dx = (maxX - minX) / 20

    y = minY + dy
    x = maxX + dx

    m = max(list(set(polynomial_orders)))

    for k,c in colors.items():
        if k <= m :
            p = Rectangle(xy=(x,y), width=dx, height=dy, fill=True, facecolor=c)
            sp.add_patch(p)
            sp.text(x + dx + (dx/2), y + (dy/4), str(k))
            y += dy
        else:
            break

    sp.text(x, y + (dy/2), str('Orders'))

    sp.set_title("Mesh")
    sp.set_aspect("equal")
    sp.autoscale_view()
    return sp.figure

def plot_mesh_mpl_simple(nodes, elements, orders=None, colors=None, axes=None,
        plot_nodes=True):
    """
    This plots the mesh using simple mpl plot commands.
    """
    if axes is None:
        from pylab import gca
        axes = gca()

    #if orders is None:
    #    orders = [1] * len(elements)

    if colors is None:
        colors = {0: '#000000', 1: '#000684', 2: '#3250fc',
            3: '#36c4ee', 4: '#04eabc', 5: '#62ff2a', 6: '#fdff07',
            7: '#ffa044', 8: '#ff1111', 9: '#b02c2c', 10: '#820f97'}

    # check that if orders and elements match (if orders are passed in)
    if orders is not None:
        assert len(elements) == len(orders)

    # join nodes with lines:
    for i, e in enumerate(elements):
        x_avg = 0
        y_avg = 0
        for k in range(len(e)):
            n1 = e[k-1]
            n2 = e[k]
            x1, y1 = nodes[n1]
            x2, y2 = nodes[n2]
            x_avg += x2
            y_avg += y2
            axes.plot([x1, x2], [y1, y2], "-",
                    color=(0, 0, 150/255.), lw=2)
        x_avg /= len(e)
        y_avg /= len(e)
        if orders:
            axes.text(x_avg, y_avg, str(orders[i]))

    if plot_nodes:
        # plot nodes:
        if isinstance(nodes, dict):
            nodes = nodes.values()
        for n in nodes:
            x = n[0]
            y = n[1]
            axes.plot([x], [y], 's', color=(150/255., 0, 0))
    return axes.figure


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

        Example:

        >>> 1 + 1
        2
        >>> 1 + 2
        3
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
                engine = mlab.get_engine()
                image = engine.current_scene
                image.scene.background = (1.0, 1.0, 1.0)
                image.scene.foreground = (0.0, 0.0, 0.0)
                mlab.colorbar(orientation="vertical")
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

    def show(self, mesh, show=True, lib="glut", notebook=False, space=None,
            filename="a.png", **options):
        if lib == "glut":
            from _hermes2d import MeshView
            m = MeshView(self._name, self._x, self._y, self._w, self._h)
            m.show(mesh)
            m.wait()
        elif lib == "mpl":
            p = plot_hermes_mesh_mpl(mesh, space=space, **options)
            if show:
                if notebook:
                    p.savefig(filename)
                else:
                    p.show()
                    import pylab
                    pylab.show()
            return p
        else:
            raise NotImplementedError("Unknown library '%s'" % lib)
