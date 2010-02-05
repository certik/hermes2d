from pylab import plot, show, legend, axes, figure, title
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

from numpy.linalg import solve
from numpy import arange, array, dot

def run(sys):
    A = sys.get_matrix().todense()
    rhs = sys.get_rhs()
    print A
    print "-"*70
    print "rhs", rhs
    x = array([
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,1,1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,1,1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,

        1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
        1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
        1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
        1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
        ])
    #x = array([
    #    1, 1, 1, 1,
    #    0, 0, 0, 0,
    #    0, 0, 0, 0,
    #    1, 0, 0, 0,])
    print x
    new_rhs = array(dot(A, x))[0]
    print new_rhs
    r = new_rhs-rhs
    print "difference"
    print r

    print "solution:"
    print solve(A, rhs)
    print solve(A, new_rhs)
    stop
    #show()

def x_reduce(x, prec=6):
    """
    Removes all duplicates and sorts the list.
    """
    x = [round(_x, prec) for _x in x]
    x = array(list(set(x)))
    x.sort()
    return x

def xyx2zformat(x, y, z):
    """
    Converts the xyz format to the mesh grid suitable for visualization.

    The x, y and z arrays have the same length and the (x[i], y[i], z[i]) tuple
    (for i in range(len(x))) contain the x, y, z coordinates of the points. The
    x, y coordinates are expected to form a 2D grid.

    xyx2zformat returns X, Y and Z, where X, Y have the same format as if
    returned by meshgrid() and Z has the same format as if returned by
    Z=sqrt(X**2+Y**2)  (just an example, any function would work).

    Example:

    >>> x, y, z = xyx2zformat([1, 1, 2, 2], [0, 1, 0, 1], [6, 7, 8, 9])
    >>> x
    array([[ 1.,  2.],
           [ 1.,  2.]])
    >>> y
    array([[ 0.,  0.],
           [ 1.,  1.]])
    >>> z
    array([[ 6.,  8.],
           [ 7.,  9.]])

    """
    prec = 6
    d = {}
    for _x, _y, _z in zip(x, y, z):
        _x = round(_x, prec)
        _y = round(_y, prec)
        d[_x] = d.get(_x, {})
        d[_x][_y] = _z
    x = x_reduce(x, prec=prec)
    y = x_reduce(y, prec=prec)
    X, Y = np.meshgrid(x, y)
    Z = X.copy()
    for i in range(len(x)):
        for j in range(len(y)):
            Z[i, j] = d[X[i, j]][Y[i, j]]
    return X, Y, Z


def plotxy(x, y, v, i):
    #if i > 4:
    #    return
    #figure()
    #title("iter=%d" % i)
    #plot(x, y, "x", label="iter=%d" % i)
    #legend()
    #axes().set_aspect('equal')
    fig = figure()
    ax = Axes3D(fig)
    X, Y, Z = xyx2zformat(x, y, v)
    #print X
    #print Y
    #print Z
    #stop
    #X, Y = np.meshgrid(X, Y)
    #R = np.sqrt(X**2 + Y**2)
    #Z = np.sin(R)
    if len(x_reduce(Z.flat)) == 1:
        # hack to make mpl plot this function:
        Z[0, 0] = 1.1
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet)
    ax.set_aspect('equal')
