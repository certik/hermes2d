Benchmarks
==========

This section contains the description of selected `benchmarks 
<http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/benchmarks>`_.
Contrary to regular examples, benchmarks typically do not have 
a strong physical or engineering motivation but they come with 
a known exact solution and thus they are a great resource for 
comparisons of various methods and adaptivity algorithms.

L-shape
-------

Complete information to this example can be found in the corresponding 
`main.cpp <http://hpfem.org/git/gitweb.cgi/hermes2d.git/blob/HEAD:/benchmarks/lshape/main.cpp>`_ file.
This is a standard adaptivity benchmark whose exact solution is smooth but
contains singular gradient in a re-entrant corner. Solved is the Laplace equation 

.. math::
    :label: lshape

       -\Delta u = 0.

Domain of interest:

.. image:: img/lshape/domain.png
   :align: center
   :width: 470
   :height: 470
   :alt: Computational domain.

Exact solution:

.. math::
    :label: lshape-exact

    u(x, y) = r^{2/3}\sin(2a/3 + \pi/3)

where $r(x,y) = \sqrt{x^2 + y^2}$ and $a(x,y) = \mbox{atan}(x/y)$. Let us show the 
computer code for the exact solution and the weak forms:

::

    // exact solution
    static double fn(double x, double y)
    {
      double r = sqrt(x*x + y*y);
      double a = atan2(x, y);
      return pow(r, 2.0/3.0) * sin(2.0*a/3.0 + M_PI/3);
    }

    static double fndd(double x, double y, double& dx, double& dy)
    {
      double t1 = 2.0/3.0*atan2(x, y) + M_PI/3;
      double t2 = pow(x*x + y*y, 1.0/3.0);
      double t3 = x*x * ((y*y)/(x*x) + 1);
      dx = 2.0/3.0*x*sin(t1)/(t2*t2) + 2.0/3.0*y*t2*cos(t1)/t3;
      dy = 2.0/3.0*y*sin(t1)/(t2*t2) - 2.0/3.0*x*t2*cos(t1)/t3;
      return fn(x, y);
    }

    // boundary condition types
    int bc_types(int marker)
    {
      return BC_ESSENTIAL;
    }

    // bilinear form corresponding to the Laplace equation
    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

Sample results of this computation are shown below.

Solution:

.. image:: img/lshape/sol_3d_view.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution.

Final mesh (h-FEM with linear elements):

.. image:: img/lshape/mesh-h1.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh (h-FEM with linear elements).

Final mesh (h-FEM with quadratic elements):

.. image:: img/lshape/mesh-h2.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. image:: img/lshape/mesh-hp.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh (hp-FEM).

DOF convergence graphs:

.. image:: img/lshape/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: img/lshape/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.




Layer
-----

To be added soon.

Line-sing
---------

To be added soon.

Kellogg
-------

To be added soon.

Smooth-iso
----------

To be added soon.

Smooth-aniso-x
--------------

To be added soon.

Smooth-aniso-y
--------------

To be added soon.

Singpert-aniso
--------------

To be added soon.

Bessel
------

To be added soon.

Screen
------

To be added soon.




























