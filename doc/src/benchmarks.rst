Benchmarks with Known Exact Solution
====================================

This section contains the description of selected `benchmarks 
<http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/benchmarks>`_.
Contrary to regular examples, benchmarks typically do not have 
a significant physical or engineering motivation but they come with 
a known exact solution and thus they are a great resource for 
comparisons of various methods and adaptivity algorithms.

Smooth-iso (Elliptic)
---------------------

**Git reference:** Benchmark `smooth-iso <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/benchmarks/smooth-iso>`_.

We show that it is a very bad idea to approximate smooth solutions using low-order 
elements.

Equation solved: Poisson equation 

.. math::
    :label: smooth-iso

       -\Delta u = f.

Domain of interest: Square $(0, \pi)^2$.

Right-hand side:

.. math::
    :label: smooth-iso-rhs
 
    f(x, y) = 2\sin(x)\sin(y).

Boundary conditions: Zero Dirichlet. 

Exact solution:

.. math::
    :label: smooth-iso-exact

    u(x, y) = \sin(x)\sin(y).

Code for the exact solution and the weak forms:

::

    // exact solution
    static double fn(double x, double y)
    {
      return sin(x)*sin(y);
    }

    static double fndd(double x, double y, double& dx, double& dy)
    {
      dx = cos(x)*sin(y);
      dy = sin(x)*cos(y);
      return fn(x, y);
    }

    // boundary condition types
    int bc_types(int marker)
    {
      return BC_ESSENTIAL;
    }

    // function values for Dirichlet boundary conditions
    scalar bc_values(int marker, double x, double y)
    {
      return fn(x, y);
    }

    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    template<typename Real>
    Real rhs(Real x, Real y)
    {
      return 2 * sin(x) * sin(y);
    }

    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_F_v<Real, Scalar>(n, wt, rhs, v, e);
    }

Solution:

.. image:: img/smooth-iso/sol_3d_view.png
   :align: center
   :width: 500
   :height: 300
   :alt: Solution.

Below we show meshes obtained using various types of adaptivity. 
Note the tremendous differences in their performance. The meshes do not correspond to 
the same level of accuracy since the low-order methods could not achieve the same error 
as hp-FEM. Therefore, compare not only the number of DOF but also the error level. 
Convergence graphs for all cases are shown at the end of this section.

Final mesh (h-FEM, p=1): 27469 DOF, error 0.39173795799476 %

.. image:: img/smooth-iso/mesh-h1.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh

Final mesh (h-FEM, p=2): 39185 DOF, error 0.0022127484879974 %

.. image:: img/smooth-iso/mesh-h2.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh

Final mesh (hp-FEM): 49 DOF, error 4.2775412425017e-05 %

.. image:: img/smooth-iso/mesh-hp.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh

DOF convergence graphs:

.. image:: img/smooth-iso/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: img/smooth-iso/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

Smooth-aniso-x (Elliptic)
-------------------------

**Git reference:** Benchmark `smooth-aniso-x <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/benchmarks/smooth-aniso-x>`_.

We show that one should use (spatially as well as polynomially) anisotropic refinements for solutions 
containing anisotropy. 

Equation solved: Poisson equation 

.. math::
    :label: sin

       -\Delta u = f.

Domain of interest: Square $(0, \pi)^2$.

Right-hand side:

.. math::
    :label: sin-rhs
 
    f(x, y) = \sin(x).

Boundary conditions: Zero Dirichlet on the left and right edges, zero Neumann on the rest of the boundary.

Exact solution:

.. math::
    :label: sin-exact

    u(x, y) = \sin(x).

Solution:

.. image:: img/smooth-aniso-x/sol_3d_view.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution.

Below we show meshes obtained using various types of adaptivity. 
Note the tremendous differences in their performance. The meshes do not correspond to 
the same level of accuracy since the low-order methods could not achieve the same error 
as hp-FEM. Therefore, compare not only the number of DOF but also the error level. 
Convergence graphs for all cases are shown at the end of this section.

Final mesh (h-FEM, p=1, isotropic refinements): 41033 DOF, error 0.22875054074711 %

.. image:: img/smooth-aniso-x/mesh-h1-iso.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh

Final mesh (h-FEM, p=1, anisotropic refinements): 39594 DOF, error 0.0039444224349215 %

.. image:: img/smooth-aniso-x/mesh-h1-aniso.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh

Final mesh (h-FEM, p=2, isotropic refinements): 54627 DOF, error 0.0017755772528929 %

.. image:: img/smooth-aniso-x/mesh-h2-iso.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh

Final mesh (h-FEM, p=2, anisotropic refinements): 3141 DOF, error 9.3084842840514e-05 %

.. image:: img/smooth-aniso-x/mesh-h2-aniso.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh

Final mesh (hp-FEM, isotropic refinements): 63 DOF, error = 3.6797337289125e-05 %

.. image:: img/smooth-aniso-x/mesh-hp-iso.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh

Final mesh (hp-FEM, anisotropic refinements): 14 DOF, error 3.6797337292196e-05 %, The 
color pattern means that the polynomial degrees are one and eight in the vertical and 
horizontal directions, respectively.

.. image:: img/smooth-aniso-x/mesh-hp-aniso.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh

DOF convergence graphs:

.. image:: img/smooth-aniso-x/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: img/smooth-aniso-x/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.


Smooth-aniso-y (Elliptic)
-------------------------

**Git reference:** Benchmark `smooth-aniso-y <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/benchmarks/smooth-aniso-y>`_.

This example is very similar to the previous one, except now the solution is 
constant in the x-direction. It is good to have both to be able to check that 
anisotropic refinements work correctly. 

L-Shape (Elliptic)
------------------

**Git reference:** Benchmark `lshape <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/benchmarks/lshape>`_.

This is a standard adaptivity benchmark whose exact solution is smooth but
contains singular gradient in a re-entrant corner. 

Equation solved: Laplace equation 

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

where $r(x,y) = \sqrt{x^2 + y^2}$ and $a(x,y) = \mbox{atan}(x/y)$. 

Code for the exact solution and the weak forms:

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

Internal Layer (Elliptic)
-------------------------

**Git reference:** Benchmark `layer <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/benchmarks/layer>`_.

This example has a smooth solution that exhibits a steep internal layer inside the domain.

Equation solved: Poisson equation 

.. math::
    :label: layer

       -\Delta u = f.

Domain of interest: Unit square $(0, 1)^2$.

Right-hand side:

.. math::
    :label: layer-rhs
 
    f(x, y) = \frac{27}{2} (2y + 0.5)^2 (\pi - 3t) \frac{S^3}{u^2 t_2} +
    \frac{27}{2} (2x - 2.5)^2 (\pi - 3t) \frac{S^3}{u^2 t_2}
    - \frac{9}{4} (2y + 0.5)^2 \frac{S}{u t^3} -
    \frac{9}{4} (2x - 2.5)^2 \frac{S}{u t^3} +
    18 \frac{S}{ut}.

Exact solution:

.. math::
    :label: layer-exact

    u(x, y) = \mbox{atan}\left(S \sqrt{(x-1.25)^2 + (y+0.25)^2} - \pi/3\right).

where $S$ is a parameter (slope of the layer). With larger $S$, this problem 
becomes difficult for adaptive algorithms, and at the same time the advantage of 
adaptive $hp$-FEM over adaptive low-order FEM becomes more significant. We will 
use $S = 60$ in the following.

Code for the exact solution and the weak forms:

::

    // exact solution
    static double fn(double x, double y)
    {
      return atan(SLOPE * (sqrt(sqr(x-1.25) + sqr(y+0.25)) - M_PI/3));
    }
    
    static double fndd(double x, double y, double& dx, double& dy)
    {
      double t = sqrt(sqr(x-1.25) + sqr(y+0.25));
      double u = t * (sqr(SLOPE) * sqr(t - M_PI/3) + 1);
      dx = SLOPE * (x-1.25) / u;
      dy = SLOPE * (y+0.25) / u;
      return fn(x, y);
    }
    
    // boundary condition types
    int bc_types(int marker)
    {
      return BC_ESSENTIAL;
    }
    
    // Dirichlet boundary condition values
    scalar bc_values(int marker, double x, double y)
    {
      return fn(x, y);
    }
    
    // bilinear form for the Poisson equation
    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }
    
    template<typename Real>
    Real rhs(Real x, Real y)
    {
      Real t2 = sqr(y + 0.25) + sqr(x - 1.25);
      Real t = sqrt(t2);
      Real u = (sqr(M_PI - 3.0*t)*sqr(SLOPE) + 9.0);
      return 27.0/2.0 * sqr(2.0*y + 0.5) * (M_PI - 3.0*t) * pow(SLOPE,3.0) / (sqr(u) * t2) +
             27.0/2.0 * sqr(2.0*x - 2.5) * (M_PI - 3.0*t) * pow(SLOPE,3.0) / (sqr(u) * t2) -
             9.0/4.0 * sqr(2.0*y + 0.5) * SLOPE / (u * pow(t,3.0)) -
             9.0/4.0 * sqr(2.0*x - 2.5) * SLOPE / (u * pow(t,3.0)) +
             18.0 * SLOPE / (u * t);
    }
     
    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return -int_F_v<Real, Scalar>(n, wt, rhs, v, e);
    }

Solution:

.. image:: img/layer/sol_3d_view.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution.

Final mesh (h-FEM with linear elements):

.. image:: img/layer/mesh-h1.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh (h-FEM with linear elements).

Final mesh (h-FEM with quadratic elements):

.. image:: img/layer/mesh-h2.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. image:: img/layer/mesh-hp.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh (hp-FEM).

DOF convergence graphs:

.. image:: img/layer/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: img/layer/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

Boundary Layer (Elliptic)
-------------------------

**Git reference:** Benchmark `layer-2 <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/benchmarks/layer-2>`_.

This example is a singularly perturbed problem with known exact solution that exhibits a thin boundary layer, that 
the reader can use to perform various experiments with adaptivity for problems with boundary layers. The sample 
numerical results presented below imply that:

* one should always use anisotropically refined meshes for problems with boundary layers,
* hp-FEM is vastly superior to h-FEM with linear and quadratic elements, 
* one should use not only spatially anisotropic elements, but also polynomial anisotropy (different polynomial orders in each direction) for problems in boundary layers. 

Equation solved: Poisson equation 

.. math::
    :label: layer-2

       -\Delta u + K^2 u = f.

Domain of interest: Square $(-1, 1)^2$.

Exact solution: 

.. math::

    u(x,y) = \hat u(x) \hat u(y)

where $\hat u$ is the exact solution of the 1D singularly-perturbed problem

.. math::

    -u'' + K^2 u = K^2

in $(-1,1)$ with zero Dirichlet boundary conditions. This solution has the form 

.. math::

    \hat u (x) = 1 - [exp(Kx) + exp(-Kx)] / [exp(K) + exp(-K)];

Right-hand side: Calculated by inserting the exact solution into the equation. Here
is the code snippet with both the exact solution and the right-hand side:

::

    // solution to the 1D problem -u'' + K*K*u = K*K in (-1,1) with zero Dirichlet BC
    double uhat(double x) {
      return 1. - (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
    }
    double duhat_dx(double x) {
      return -K * (exp(K*x) - exp(-K*x)) / (exp(K) + exp(-K));
    }
    double dduhat_dxx(double x) {
      return -K*K * (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
    }

    // exact solution u(x,y) to the 2D problem is defined as the
    // Cartesian product of the 1D solutions
    static double sol_exact(double x, double y, double& dx, double& dy)
    {
      dx = duhat_dx(x) * uhat(y);
      dy = uhat(x) * duhat_dx(y);
      return uhat(x) * uhat(y);
    }

    // right-hand side
    double rhs(double x, double y) {
      return -(dduhat_dxx(x)*uhat(y) + uhat(x)*dduhat_dxx(y)) + K*K*uhat(x)*uhat(y);
    }

The weak forms are very simple and they are defined as follows. The only thing worth mentioning 
here is that we integrate the non-polynomial right-hand side with a very hign order for accuracy:

::

    // weak forms
    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) + K*K * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_F_v<Real, Scalar>(n, wt, rhs, v, e);;
    }

    // integration order for linear_form_0
    Ord linear_form_ord(int n, double *wt, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return 24;
    }

The numerical results follow:

Solution:

.. image:: img/layer-2/solution.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution.

Below we present a series of convergence comparisons. Note that the error plotted
is the true approximate error calculated wrt. the exact solution given above.

Isotropic refinements
~~~~~~~~~~~~~~~~~~~~~

Let us first compare the performance of h-FEM (p=1), h-FEM (p=2) and hp-FEM with **isotropic** refinements:

Final mesh (h-FEM, p=1, isotropic refinements):

.. image:: img/layer-2/mesh_h1_iso.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh.

Final mesh (h-FEM, p=2, isotropic refinements):

.. image:: img/layer-2/mesh_h2_iso.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh.

Final mesh (hp-FEM, isotropic refinements):

.. image:: img/layer-2/mesh_hp_iso.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh.

DOF convergence graphs:

.. image:: img/layer-2/conv_compar_dof_iso.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU convergence graphs:

.. image:: img/layer-2/conv_compar_cpu_iso.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

Anisotropic refinements
~~~~~~~~~~~~~~~~~~~~~~~

Next we compare the performance of h-FEM (p=1), h-FEM (p=2) and hp-FEM with **anisotropic** refinements:

Final mesh (h-FEM, p=1, anisotropic refinements):

.. image:: img/layer-2/mesh_h1_aniso.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh.

Final mesh (h-FEM, p=2, anisotropic refinements):

.. image:: img/layer-2/mesh_h2_aniso.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh.

Final mesh (hp-FEM, anisotropic refinements):

.. image:: img/layer-2/mesh_hp_aniso.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh.

DOF convergence graphs:

.. image:: img/layer-2/conv_compar_dof_aniso.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU convergence graphs:

.. image:: img/layer-2/conv_compar_cpu_aniso.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

h-FEM (p=1): comparison of isotropic and anisotropic refinements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DOF convergence graphs:

.. image:: img/layer-2/conv_compar_dof_h1.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU convergence graphs:

.. image:: img/layer-2/conv_compar_cpu_h1.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

h-FEM (p=2): comparison of isotropic and anisotropic refinements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DOF convergence graphs:

.. image:: img/layer-2/conv_compar_dof_h2.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU convergence graphs:

.. image:: img/layer-2/conv_compar_cpu_h2.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

hp-FEM: comparison of isotropic and anisotropic refinements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the hp-FEM one has two kinds of anisotropy -- spatial and polynomial. In the following,
"iso" means isotropy both in h and p, "aniso h" means anisotropy in h only, and 
"aniso hp" means anisotropy in both h and p. 

DOF convergence graphs (hp-FEM):

.. image:: img/layer-2/conv_compar_dof_hp.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU convergence graphs (hp-FEM):

.. image:: img/layer-2/conv_compar_cpu_hp.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

The reader can see that enabling polynomially anisotropic refinements in the hp-FEM is 
equally important as allowing spatially anisotropic ones. 

Line singularity (Elliptic)
---------------------------

**Git reference:** Benchmark `line-singularity <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/benchmarks/line-singularity>`_.

The is another example with anisotropic solution that is suitable for testing 
anisotropic element refinements.

Equation solved: Poisson equation 

.. math::
    :label: line-sing

       -\Delta u = f.

Domain of interest: Square $(-1, 1)^2$.

Boundary conditions: Zero Neumann on left edge, Dirichlet given by the 
exact solution on the rest of the boundary.

Exact solution: 

.. math::

    u(x,y) = \cos(Ky)\ \ \ \mbox{for}\ x \le 0,\\
    u(x,y) = \cos(Ky) + x^{\alpha}\ \ \ \mbox{for}\ x > 0,

where $K$ and $\alpha$ are real constants. 

Right-hand side: Obtained by inserting the exact solution into the equation.
The corresponding code snippet is shown below:

::

    scalar rhs(scalar x, scalar y)
    {
      if (x < 0) return fn(x, y)*K*K;
      else return fn(x, y)*K*K-ALPHA*(ALPHA-1)*pow(x, ALPHA - 2.) - K*K*pow(x, ALPHA);
    }

Solution for $K = \pi/2$ and $\alpha = 2.01$:

.. image:: img/line-singularity/solution.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution.

Comparison of h-FEM (p=1), h-FEM (p=2) and hp-FEM with anisotropic refinements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Final mesh (h-FEM, p=1, anisotropic refinements):

.. image:: img/line-singularity/mesh_h1_aniso.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (h-FEM, p=2, anisotropic refinements):

.. image:: img/line-singularity/mesh_h2_aniso.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (hp-FEM, h-anisotropic refinements):

.. image:: img/line-singularity/mesh_hp_anisoh.png
   :align: center
   :width: 450
   :alt: Final mesh.

DOF convergence graphs:

.. image:: img/line-singularity/conv_dof_aniso.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU convergence graphs:

.. image:: img/line-singularity/conv_cpu_aniso.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

hp-FEM with iso, h-aniso and hp-aniso refinements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Final mesh (hp-FEM, isotropic refinements):

.. image:: img/line-singularity/mesh_hp_iso.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (hp-FEM, h-anisotropic refinements):

.. image:: img/line-singularity/mesh_hp_anisoh.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (hp-FEM, hp-anisotropic refinements):

.. image:: img/line-singularity/mesh_hp_aniso.png
   :align: center
   :width: 450
   :alt: Final mesh.

The following convergence comparisons show that, surprizingly, 
the h-anisotropic version is the most efficient one. The reason is
that both the hp-FEM with isotropic refinements and hp-FEM with hp-anisotropic 
refinements hit the maximum polynomial degree (p=9) early in the mesh 
adaptation process (while hp-FEM with h-anisotropic refinements does not). 
We expect that this problem should go away after we increase the 
maximum polynomial degree in Hermes - this is on the TODO list now.  

DOF convergence graphs:

.. image:: img/line-singularity/conv_dof_hp.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU convergence graphs:

.. image:: img/line-singularity/conv_cpu_hp.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.



Kellogg (Elliptic)
------------------

**Git reference:** Benchmark `kellogg <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/benchmarks/kellogg>`_.

The solution to this elliptic problems contains a severe singularity that poses a challenge to 
adaptive methods. 

Equation solved:

.. math::

       -\nabla \cdot (a(x,y) \nabla u) = 0,

where the parameter $a$ is piecewise-constant, $a(x,y) = R$ in the first and third quadrants and $a(x,y) = 1$ 
in the remaining two quadrants. 

Domain of interest: Square $(-1, 1)^2$.

Right-hand side: $f(x,y) = 0$.

Boundary conditions: Dirichlet given by exact solution. 

Exact solution: Quite complicated, see the code below.

::

    // problem constants
    const double R = 161.4476387975881;      // Equation parameter.
    const double TAU = 0.1;                  // Equation parameter.
    const double RHO = M_PI/4.;              // Equation parameter
    const double SIGMA = -14.92256510455152; // Equation parameter

    // exact solution
    static double fn(double x, double y)
    {
      double theta = atan2(y,x);
      if (theta < 0) theta = theta + 2.*M_PI;
      double r = sqrt(x*x + y*y);

      double mu;
      if (theta <= M_PI/2.) {
        mu = cos((M_PI/2. - SIGMA)*TAU) * cos((theta - M_PI/2. + RHO)*TAU);
      }
      else {
        if (theta <= M_PI) {
          mu = cos(RHO*TAU) * cos((theta - M_PI + SIGMA)*TAU);
        }
        else {
          if (theta <= 3.*M_PI/2.) {
            mu = cos(SIGMA*TAU) * cos((theta - M_PI - RHO)*TAU);
          }
          else {
            mu = cos((M_PI/2. - RHO)*TAU) * cos((theta - 3.*M_PI/2. - SIGMA)*TAU);
          }
        }
      }

      return pow(r, TAU) * mu;
    }

The weak forms are as follows:

::

    // Weak forms
    template<typename Real, typename Scalar>
    Scalar bilinear_form_I_III(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return R*int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    template<typename Real, typename Scalar>
    Scalar bilinear_form_II_IV(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return 1.*int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }


Solution:

.. image:: img/kellogg/solution.png
   :align: center
   :width: 600
   :alt: Solution.

Final mesh (h-FEM with linear elements):

.. image:: img/kellogg/mesh-h1.png
   :align: center
   :width: 600
   :alt: Mesh.

Final mesh (h-FEM with quadratic elements):

.. image:: img/kellogg/mesh-h2.png
   :align: center
   :width: 600
   :alt: Mesh.

Final mesh (hp-FEM):

.. image:: img/kellogg/mesh-hp.png
   :align: center
   :width: 600
   :alt: Mesh.

DOF convergence graphs:

.. image:: img/kellogg/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: img/kellogg/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

Bessel (Maxwell's Equations)
----------------------------

**Git reference:** Benchmark `bessel <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/benchmarks/bessel>`_.

This example solves time-harmonic Maxwell's equations in an L-shaped domain and it 
describes the diffraction of an electromagnetic wave from a re-entrant corner. It comes with an 
exact solution that contains singularity.

Equation solved: Time-harmonic Maxwell's equations

.. math::
    :label: bessel

    \frac{1}{\mu_r} \nabla \times \nabla \times E - \kappa^2 \epsilon_r E = \Phi.

Domain of interest is the square $(-10, 10)^2$ missing the quarter lying in the 
fourth quadrant. It is filled with air:

.. image:: img/bessel/domain.png
   :align: center
   :width: 490
   :height: 490
   :alt: Computational domain.

Boundary conditions: Combined essential and natural, see the 
`main.cpp <http://hpfem.org/git/gitweb.cgi/hermes2d.git/blob/HEAD:/benchmarks/bessel/main.cpp>`_ file.

Exact solution:

.. math::
    :label: bessel-exact

    E(x, y) = \nabla \times J_{\alpha} (r) \cos(\alpha \theta)

where $J_{\alpha}$ is the Bessel function of the first kind, 
$(r, \theta)$ the polar coordinates and $\alpha = 2/3$. In 
computer code, this reads:

::

    void exact_sol(double x, double y, scalar& e0, scalar& e1)
    {
      double t1 = x*x;
      double t2 = y*y;
      double t4 = sqrt(t1+t2);
      double t5 = jv(-1.0/3.0,t4);
      double t6 = 1/t4;
      double t7 = jv(2.0/3.0,t4);
      double t11 = (t5-2.0/3.0*t6*t7)*t6;
      double t12 = atan2(y,x);
      if (t12 < 0) t12 += 2.0*M_PI;
      double t13 = 2.0/3.0*t12;
      double t14 = cos(t13);
      double t17 = sin(t13);
      double t18 = t7*t17;
      double t20 = 1/t1;
      double t23 = 1/(1.0+t2*t20);
      e0 = t11*y*t14-2.0/3.0*t18/x*t23;
      e1 = -t11*x*t14-2.0/3.0*t18*y*t20*t23;
    }  

Here jv() is the Bessel function $\bfJ_{\alpha}$. For its source code see the 
`forms.cpp <http://hpfem.org/git/gitweb.cgi/hermes2d.git/blob/HEAD:/benchmarks/bessel/forms.cpp>`_ file.

Code for the weak forms:

::

    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
    return 1.0/mu_r * int_curl_e_curl_f<Real, Scalar>(n, wt, u, v) -
           sqr(kappa) * int_e_f<Real, Scalar>(n, wt, u, v);
    }
   
    template<typename Real, typename Scalar>
    Scalar bilinear_form_surf(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      cplx ii = cplx(0.0, 1.0);
      return ii * (-kappa) * int_e_tau_f_tau<Real, Scalar>(n, wt, u, v, e);
    }
   
    scalar linear_form_surf(int n, double *wt, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      scalar result = 0;
      for (int i = 0; i < n; i++)
      {
        double r = sqrt(e->x[i] * e->x[i] + e->y[i] * e->y[i]);
        double theta = atan2(e->y[i], e->x[i]);
        if (theta < 0) theta += 2.0*M_PI;
        double j13    = jv(-1.0/3.0, r),    j23    = jv(+2.0/3.0, r);
        double cost   = cos(theta),         sint   = sin(theta);
        double cos23t = cos(2.0/3.0*theta), sin23t = sin(2.0/3.0*theta);
   
        double Etau = e->tx[i] * (cos23t*sint*j13 - 2.0/(3.0*r)*j23*(cos23t*sint + sin23t*cost)) +
                      e->ty[i] * (-cos23t*cost*j13 + 2.0/(3.0*r)*j23*(cos23t*cost - sin23t*sint));
  
        result += wt[i] * cplx(cos23t*j23, -Etau) * ((v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
      }
      return result;
    }
    // maximal polynomial order to integrate surface linear form
    Ord linear_form_surf_ord(int n, double *wt, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {  return Ord(v->val[0].get_max_order());  }

Solution:

.. image:: img/bessel/solution.png
   :align: center
   :width: 500
   :height: 420
   :alt: Solution.

Final mesh (h-FEM with linear elements):

.. image:: img/bessel/mesh-h1.png
   :align: center
   :width: 460
   :height: 390
   :alt: Final mesh (h-FEM with linear elements).

Note that the polynomial order indicated corresponds to the tangential components 
of approximation on element interfaces, not to polynomial degrees inside the elements
(those are one higher).

Final mesh (h-FEM with quadratic elements):

.. image:: img/bessel/mesh-h2.png
   :align: center
   :width: 460
   :height: 390
   :alt: Final mesh (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. image:: img/bessel/mesh-hp.png
   :align: center
   :width: 460
   :height: 390
   :alt: Final mesh (hp-FEM).

DOF convergence graphs:

.. image:: img/bessel/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: img/bessel/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

Screen (Maxwell's Equations)
----------------------------

**Git reference:** Benchmark `screen <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/benchmarks/screen>`_.

This example solves time-harmonic Maxwell's equations. It describes an electromagnetic wave that 
hits a thin screen under the angle of 45 degrees, causing a singularity at the tip of the screen.
The strength of the singularity makes this example rather difficult. 

Equation solved: Time-harmonic Maxwell's equations

.. math::
    :label: screen

    \frac{1}{\mu_r} \nabla \times \nabla \times E - \kappa^2 \epsilon_r E = \Phi.

Domain of interest is the square $(-1,1)^2$ missing the edge that connects the center with 
the midpoint of the left side. It is filled with air:

.. image:: img/screen/domain.png
   :align: center
   :width: 490
   :height: 490
   :alt: Computational domain.

Boundary conditions: Tangential component of solution taken from known exact solution 
(essential BC). See the 
`main.cpp <http://hpfem.org/git/gitweb.cgi/hermes2d.git/blob/HEAD:/benchmarks/screen/main.cpp>`_ file.

Exact solution: This is rather complicated in this case - please look into the 
corresponding file 
`exact_sol.cpp <http://hpfem.org/git/gitweb.cgi/hermes2d.git/blob/HEAD:/benchmarks/screen/exact_sol.cpp>`_.

Code for the weak forms:

::

    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_curl_e_curl_f<Real, Scalar>(n, wt, u, v) - int_e_f<Real, Scalar>(n, wt, u, v);
    }

Solution (real part of $E_1$):

.. image:: img/screen/sol1.png
   :align: center
   :width: 510
   :height: 470
   :alt: Solution.

Solution (real part of $E_2$):

.. image:: img/screen/sol2.png
   :align: center
   :width: 510
   :height: 470
   :alt: Solution.

Solution (imaginary part of $E_1$):

.. image:: img/screen/sol3.png
   :align: center
   :width: 510
   :height: 470
   :alt: Solution.

Solution (imaginary part of $E_2$):

.. image:: img/screen/sol4.png
   :align: center
   :width: 510
   :height: 470
   :alt: Solution.

Final mesh (h-FEM with linear elements):

.. image:: img/screen/mesh-h1.png
   :align: center
   :width: 460
   :height: 410
   :alt: Final mesh (h-FEM with linear elements).

Note that the polynomial order indicated corresponds to the tangential components 
of approximation on element interfaces, not to polynomial degrees inside the elements
(those are one higher).

Final mesh (h-FEM with quadratic elements):

.. image:: img/screen/mesh-h2.png
   :align: center
   :width: 460
   :height: 410
   :alt: Final mesh (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. image:: img/screen/mesh-hp.png
   :align: center
   :width: 460
   :height: 410
   :alt: Final mesh (hp-FEM).

DOF convergence graphs:

.. image:: img/screen/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: img/screen/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.































