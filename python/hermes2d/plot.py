from hermes2d import Linearizer

def sln2png(sln, filename):
    """
    Creates a nice png image of the Solution sln.
    """
    plot_sln_mayavi(sln)
    from enthought.mayavi.mlab import savefig
    savefig(filename)

def plot_sln_mpl(sln, method="default"):
    """
    Plots the Solution() instance sln using Linearizer() and matplotlib.

    method = "default" ... creates a plot using triangles (the triangles are
                not interpolated, so sometimes one can see small defects)
    method = "contour" ... takes the vertices from linearizer and interpolates
                them using contour and contourf (it doesn't take into account
                the triangulation, so one can see the defects from the convex
                hull approximation)
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
        vals = array(vals)
        col = collections.PolyCollection(verts, linewidths=0, antialiaseds=0)
        col.set_array(vals)
        col.set_cmap(plt.cm.jet)
        fig = plt.figure()
        ax = fig.gca()
        ax.add_collection(col)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_aspect("equal")
        #plt.colorbar()
        plt.title('Solution')
    else:
        raise ValueError("Unknown method (%s)" % method)

def plot_sln_mayavi(sln):
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
    # the off screen rendering properly works only with VTK-5.2 or above:
    mlab.options.offscreen = True
    s = mlab.triangular_mesh(x, y, z, triangles, scalars=t)
    mlab.view(azimuth=90, elevation=180)
    return s

def plot_sln_pyglet(sln):
    from numpy import concatenate, array
    import pylab
    import pyglet
    from pyglet_camera import Camera
    from math import pi
    l = Linearizer()
    l.process_solution(sln)
    v = l.get_vertices()
    t = l.get_triangles()
    vertices_indices = concatenate(t)
    vertices_coordinates = concatenate(v[:, :2])
    colors = v[:, 2]
    c = pylab.get_cmap("jet")
    d = lambda x: c(x)[:3]
    colors = (colors-colors.min())/(colors.max()-colors.min())
    print "converting colors..."
    colors = [d(x) for x in colors]
    print "    done."
    colors = concatenate(colors)
    colors = array(colors*255, dtype=int)

    window = pyglet.window.Window()
    # adjust the camera so that the whole mesh is seen:
    x_length = v[:, 0].max() - v[:, 0].min()
    y_length = v[:, 1].max() - v[:, 1].min()
    max_length = max(x_length, y_length)/2.
    camera = Camera((0, 0), max_length)

    pyglet.gl.glClearColor(1, 1, 1, 1)


    @window.event
    def on_draw():
        window.clear()
        camera.update()
        camera.focus(window.width, window.height)
        pyglet.graphics.draw_indexed(len(v), pyglet.gl.GL_TRIANGLES,
                vertices_indices,
                ('v2f', vertices_coordinates),
                ('c3B', colors)
                )
        camera.hud_mode(window.width, window.height)

    pyglet.clock.schedule(lambda _: None)

    @window.event
    def on_key_press(symbol, modifiers):
        if symbol in [pyglet.window.key.Q, pyglet.window.key.ESCAPE]:
            window.close()
        elif symbol == pyglet.window.key.LEFT:
            camera.pan(camera.scale, +pi/2)
        elif symbol == pyglet.window.key.RIGHT:
            camera.pan(camera.scale, -pi/2)
        elif symbol == pyglet.window.key.UP:
            camera.pan(camera.scale, pi)
        elif symbol == pyglet.window.key.DOWN:
            camera.pan(camera.scale, 0)
        elif symbol == pyglet.window.key.COMMA:
            camera.tilt(-1)
        elif symbol == pyglet.window.key.PERIOD:
            camera.tilt(+1)
        elif symbol == pyglet.window.key.PAGEUP:
            camera.zoom(2)
        elif symbol == pyglet.window.key.PAGEDOWN:
            camera.zoom(0.5)

    @window.event
    def on_mouse_scroll(x, y, scroll_x, scroll_y):
        if scroll_y > 0:
            camera.zoom(0.5)
        else:
            camera.zoom(2)

    @window.event
    def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
        if buttons & pyglet.window.mouse.LEFT:
            camera.move(-dx, -dy)
        elif buttons & pyglet.window.mouse.RIGHT:
            scale = 1.03
            if dy > 0:
                scale = 1/scale
            camera.zoom(scale)

    print "plotting"
    pyglet.app.run()

class ScalarView(object):

    def __init__(self, x=0, y=0, w=50, h=50, name="Solution"):
        self._name = name

    def show_scale(self, *args):
        pass

    def show_mesh(self, *args):
        pass

    def wait(self):
        pass

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
        if lib == "mpl":
            plot_sln_mpl(sln, **options)
            import pylab
            if show:
                if notebook:
                    pylab.savefig(filename)
                else:
                    pylab.show()
        elif lib == "mayavi":
            plot_sln_mayavi(sln)
            from enthought.mayavi import mlab
            if show:
                if notebook:
                    mlab.savefig(filename)
                else:
                    mlab.show()
        elif lib == "pyglet":
            if show and not notebook:
                plot_sln_pyglet(sln)
            else:
                ValueError("pyglet only works with show=True and notebook=False")

        else:
            raise NotImplementedError("Unknown library '%s'" % lib)

class MeshView(object):

    def __init__(self, name="Solution", x=0, y=0, w=50, h=50):
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
        else:
            raise NotImplementedError("Unknown library '%s'" % lib)
