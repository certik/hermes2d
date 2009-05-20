from hermes2d import Linearizer

def sln2png(sln, filename):
    """
    Creates a nice png image of the Solution sln.
    """
    plot_sln_mayavi(sln)
    from enthought.mayavi.mlab import savefig
    savefig(filename)

def plot_sln_mpl(sln):
    """
    Plots the Solution() instance sln using Linearizer() and matplotlib.

    Currently only a very simple version is implemented, that takes the
    vertices from linearizer and interpolates them. More sophisticated version
    should take the triangles.
    """
    lin = Linearizer()
    lin.process_solution(sln)
    vert = lin.get_vertices()
    from numpy import min, max, linspace
    from matplotlib.mlab import griddata
    import matplotlib.pyplot as plt
    # make up data.
    #npts = int(raw_input('enter # of random points to plot:'))
    npts = 200
    x = vert[:, 0]
    y = vert[:, 1]
    z = vert[:, 2]
    # define grid.
    xi = linspace(min(x), max(x), 100)
    yi = linspace(min(y), max(y), 100)
    # grid the data.
    zi = griddata(x,y,z,xi,yi)
    # contour the gridded data, plotting dots at the nonuniform data points.
    CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
    CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
    plt.colorbar() # draw colorbar
    plt.title('Solution (%d points)' % npts)

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

class ScalarView(object):

    def __init__(self, x=0, y=0, w=50, h=50, name="Solution"):
        self._name = name

    def show_scale(self, *args):
        pass

    def show_mesh(self, *args):
        pass

    def wait(self):
        pass

    def show(self, sln, show=True, lib="mpl", notebook=False, filename="a.png"):
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
            plot_sln_mpl(sln)
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
        else:
            raise NotImplementedError("Unknown library '%s'" % lib)
