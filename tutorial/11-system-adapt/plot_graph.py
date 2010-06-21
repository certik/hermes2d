# import libraries
import pylab
#from pylab import *

# plot DOF convergence graph
pylab.axis('equal')
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
data = pylab.loadtxt("conv_dof_m.dat")
x = data[:, 0]
y = data[:, 1]
pylab.loglog(x, y, "-s", label="error (multi-mesh)")
#data = numpy.loadtxt("conv_dof_s.dat")
#x = data[:, 0]
#y = data[:, 1]
#loglog(x, y, "-s", label="error (single-mesh)")

pylab.legend()

# initialize new window
pylab.figure()

pylab.axis('equal')
pylab.title("Error convergence")
pylab.xlabel("CPU time (s)")
pylab.ylabel("Error [%]")
data = pylab.loadtxt("conv_cpu_m.dat")
x = data[:, 0]
y = data[:, 1]
pylab.loglog(x, y, "-s", label="error (multi-mesh)")
#data = numpy.loadtxt("conv_cpu_s.dat")
#x = data[:, 0]
#y = data[:, 1]
#loglog(x, y, "-s", label="error (single-mesh)")
pylab.legend()


# finalize
pylab.show()
