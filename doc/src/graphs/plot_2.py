from pylab import plot, show, legend, log
import pylab, numpy
pylab.yscale("log")
pylab.title("Error Convergence for the Bracket Problem")
pylab.xlabel("Degrees of Freedom")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_dof_m.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="hp-FEM (multi-mesh)")
data = numpy.loadtxt("conv_dof_s.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="hp-FEM (standard)")

legend()
show()
