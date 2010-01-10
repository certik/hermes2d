from pylab import plot, show, legend, log
import pylab, numpy
pylab.yscale("log")
pylab.title("Error Convergence")
pylab.xlabel("CPU Time")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_cpu_h1.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="h-FEM (p=1)")
data = numpy.loadtxt("conv_cpu_h2.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="h-FEM (p=2)")
data = numpy.loadtxt("conv_cpu_hp.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="hp-FEM")
legend()
show()
