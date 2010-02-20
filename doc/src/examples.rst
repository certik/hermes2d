Examples
========

This section contains the description of selected `examples 
<http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/examples>`_.
Its purpose is to complement rather than duplicate the information 
in the source code.

Saphir
------

More information to this example can be found in the corresponding 
`main.cpp <http://hpfem.org/git/gitweb.cgi/hermes2d.git/blob/HEAD:/examples/saphir/main.cpp>`_ file.
This is a standard nuclear engineering benchmark (IAEA number EIR-2) describing 
an external-force-driven configuration without fissile materials present, using one-group 
neutron diffusion approximation

.. math::
    :label: saphir

       -\nabla \cdot (D(x,y) \nabla \Phi) + \Sigma_a(x,y) \Phi = Q_{ext}(x,y).

The domain of interest is a 96 x 86 cm rectangle consisting of five regions:

.. image:: img/saphir/saphir.png
   :align: center
   :width: 400
   :height: 400
   :alt: Schematic picture for the saphir example.

The unknown is the neutron flux $\Phi(x, y)$. The values of the diffusion coefficient 
$D(x, y)$, absorption cross-section $\Sigma_a(x, y)$ and the source term $Q_{ext}(x,y)$
are constant in the subdomains. The source $Q_{ext} = 1$ in areas 1 and 3 and zero 
elsewhere. Boundary conditions for the flux $\Phi$ are zero everywhere. 

It is worth noticing how different material parameters are handled - we define a separate weak form 
for each material. This approach is more flexible to how material parameters were handled in 
tutorial examples 07 and 12:

::

    // Bilinear form (material 1)  
    template<typename Real, typename Scalar>
    Scalar bilinear_form_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return D_1 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) 
             + SIGMA_A_1 * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    // Bilinear form (material 2)
    template<typename Real, typename Scalar>
    Scalar bilinear_form_2(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return D_2 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) 
             + SIGMA_A_2 * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    // Bilinear form (material 3)
    template<typename Real, typename Scalar>
    Scalar bilinear_form_3(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return D_3 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) 
             + SIGMA_A_3 * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    // Bilinear form (material 4)
    template<typename Real, typename Scalar>
    Scalar bilinear_form_4(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return D_4 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) 
             + SIGMA_A_4 * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    // Bilinear form (material 5)
    template<typename Real, typename Scalar>
    Scalar bilinear_form_5(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return D_5 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) 
             + SIGMA_A_5 * int_u_v<Real, Scalar>(n, wt, u, v);
    }

The weak forms are associated with element material flags (coming from the mesh file) as follows:

::

    // initialize the weak formulation
    WeakForm wf(1);
    wf.add_biform(0, 0, bilinear_form_1, bilinear_form_ord, SYM, 1);
    wf.add_biform(0, 0, bilinear_form_2, bilinear_form_ord, SYM, 2);
    wf.add_biform(0, 0, bilinear_form_3, bilinear_form_ord, SYM, 3);
    wf.add_biform(0, 0, bilinear_form_4, bilinear_form_ord, SYM, 4);
    wf.add_biform(0, 0, bilinear_form_5, bilinear_form_ord, SYM, 5);
    wf.add_liform(0, linear_form_1, linear_form_ord, 1);
    wf.add_liform(0, linear_form_3, linear_form_ord, 3);

Sample results of this computation are shown below.

Solution:

.. image:: img/saphir/saphir-sol.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution to the saphir example.

Final mesh (h-FEM with linear elements):

.. image:: img/saphir/saphir-mesh-h1.png
   :align: center
   :width: 440
   :height: 400
   :alt: Final finite element mesh for the saphir example (h-FEM with linear elements).

Final mesh (h-FEM with quadratic elements):

.. image:: img/saphir/saphir-mesh-h2.png
   :align: center
   :width: 440
   :height: 400
   :alt: Final finite element mesh for the saphir example (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. image:: img/saphir/saphir-mesh-hp.png
   :align: center
   :width: 440
   :height: 400
   :alt: Final finite element mesh for the saphir example (hp-FEM).

DOF convergence graphs:

.. image:: img/saphir/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph for example saphir.

CPU time convergence graphs:

.. image:: img/saphir/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph for example saphir.

Iron-Water
----------

More information to this example can be found in the corresponding 
`main.cpp <http://hpfem.org/git/gitweb.cgi/hermes2d.git/blob/HEAD:/examples/iron-water/main.cpp>`_ file.
This example is very similar to the example "saphir", the main difference being that 
it reads a mesh file in the exodusii format (created by Cubit). 
This example only builds 
if you have the `ExodusII <http://sourceforge.net/projects/exodusii/>`_ and 
`NetCDF <http://www.unidata.ucar.edu/software/netcdf/>`_ libraries installed on your system and 
the variables WITH_EXODUSII, EXODUSII_ROOT and NETCDF_ROOT defined properly. 
The latter can be done, for example, in the CMake.vars file as follows:

::

    SET(WITH_EXODUSII YES)
    SET(EXODUSII_ROOT /opt/packages/exodusii)
    SET(NETCDF_ROOT   /opt/packages/netcdf)

The mesh is now loaded using the ExodusIIReader (see 
the `mesh_loader.h <http://hpfem.org/git/gitweb.cgi/hermes2d.git/blob/HEAD:/src/mesh_loader.h>`_ file):

::

    // Load the mesh
    Mesh mesh;
    ExodusIIReader mloader;
    if (!mloader.load("iron-water.e", &mesh)) error("ExodusII mesh load failed.");

The model describes an external-force-driven configuration without fissile materials present.
We will solve the one-group neutron diffusion equation

.. math::
    :label: iron-water

       -\nabla \cdot (D(x,y) \nabla \Phi) + \Sigma_a(x,y) \Phi = Q_{ext}(x,y).

The domain of interest is a 30 x 30 cm square consisting of four regions.
A uniform volumetric source is placed in water in the lower-left corner 
of the domain, surrounded with a layer of water, a layer of iron, and finally
another layer of water:

.. image:: img/iron-water/iron-water.png
   :align: center
   :width: 400
   :height: 400
   :alt: Schematic picture for the iron-water example.

The unknown is the neutron flux $\Phi(x, y)$. The values of the diffusion coefficient 
$D(x, y)$, absorption cross-section $\Sigma_a(x, y)$ and the source term $Q_{ext}(x,y)$
are constant in the subdomains. The source $Q_{ext} = 1$ in area 1 and zero 
elsewhere. The boundary conditions for this problem are zero Dirichlet (right and top edges)
and zero Neumann (bottom and left edges). Sample results of this computation are shown below.

Solution:

.. image:: img/iron-water/iron-water-sol.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution to the iron-water example.


Final mesh (h-FEM with linear elements):

.. image:: img/iron-water/iron-water-mesh-h1.png
   :align: center
   :width: 440
   :height: 400
   :alt: Final finite element mesh for the iron-water example (h-FEM with linear elements).

Final mesh (h-FEM with quadratic elements):

.. image:: img/iron-water/iron-water-mesh-h2.png
   :align: center
   :width: 440
   :height: 400
   :alt: Final finite element mesh for the iron-water example (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. image:: img/iron-water/iron-water-mesh-hp.png
   :align: center
   :width: 440
   :height: 400
   :alt: Final finite element mesh for the iron-water example (hp-FEM).

DOF convergence graphs:

.. image:: img/iron-water/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph for example iron-water.

CPU time convergence graphs:

.. image:: img/iron-water/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph for example iron-water.

Navier-Stokes
-------------

More information to this example can be found in the corresponding 
`main.cpp <http://hpfem.org/git/gitweb.cgi/hermes2d.git/blob/HEAD:/examples/ns-timedep/main.cpp>`_ file.
This model problem is concerned with the approximate solution of external
flow past a cylinder with unit diameter. The corresponding files can be found in 
example `ns-timedep <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/examples/ns-timedep>`_. The geometry looks as follows:

.. image:: img/cylinder.png
   :align: center
   :width: 600
   :height: 220
   :alt: Domain for the Navier-Stokes problem.


The motion of the fluid is described by the dimensionless incompressible
Navier-Stokes equations,

.. math::
    :label: ns1

         \dd{\bfu}{t} - \frac{1}{\rm Re} \Delta \bfu + (\bfu \cdot \nabla) \bfu + \nabla p  = 0,


.. math::
    :label: ns2

         \nabla \cdot \bfu = 0,

where $\bfu = (u_1, u_2)^T$ is the fluid velocity, $p$ is the kinematic pressure and Re
is the Reynods number. One way to solve the nonlinear system :eq:`ns1`--:eq:`ns2` is to
introduce a small time step $\tau > 0$, replace the time derivative by a backward
difference formula and linearize the convective term
$(\bfu \cdot \nabla) \bfu \approx (\bfu^{n-1} \cdot \nabla) \bfu^n$, where $\bfu^n$ is the
approximate solution on the $n$-th time level. This leads to a system of linear PDEs for the
$n$-th time level

.. math::
    :label: ns3

         \frac{\bfu^n - \bfu^{n-1}}{\tau} - \frac{1}{\rm Re} \Delta \bfu^n +     (\bfu^{n-1} \cdot \nabla) \bfu^n + \nabla p  = 0,


.. math::
    :label: ns4

         \nabla \cdot \bfu^n = 0,

Testing :eq:`ns3` by the velocity test functions $(v_1, v_2)$ and testing :eq:`ns4`
by the pressure test function $q$, we obtain the following weak formulation:

.. math::

    \int_\Omega \frac{u_1 v_1}{\tau} +   \frac{1}{\rm Re} \nabla u_1 \cdot \nabla v_1 +   (\bfu^{n-1} \cdot \nabla) u_1 v_1 - p \dd{v_1}{x} \dx   = \int_\Omega \frac{u^{n-1}_1 v_1}{\tau}


.. math::

    \int_\Omega \frac{u_2 v_2}{\tau} +   \frac{1}{\rm Re} \nabla u_2 \cdot \nabla v_2 +   (\bfu^{n-1} \cdot \nabla) u_2 v_2 - p \dd{v_2}{y} \dx   = \int_\Omega \frac{u^{n-1}_2 v_2}{\tau}


.. math::

    \int_\Omega \dd{u_1}{x} q + \dd{u_2}{y} q \dx = 0


The boundary and initial conditions for the problem are

.. math::

    \bfu(\bfx, t) = (1, 0)^T \ \ \ \ \mbox{on}\ \ \Gamma_1 \cup \Gamma_3 \cup \Gamma_4


.. math::

    \bfu(\bfx, t) = (0, 0)^T \ \ \ \ \mbox{on}\ \ \Gamma_5


.. math::

    \mbox{\it ``do-nothing"}\ \ \ \ \mbox{on}\ \ \Gamma_2


.. math::
    :label: ns:initial

     \bfu(\bfx, 0) = \bfu^0 = (0, 0)^T


In CFD, the *do-nothing* condition is a common artificial boundary condition defining
an outlet for the fluid. It means that there is no restriction on the value
of the velocity on $\Gamma_2$.

The implementation starts by defining three spaces xvel, yvel and press
for the three solution components $u_1$, $u_2$ and $p$. Using Space::set_bc_type()
we denote the Dirichlet boundary for velocity:
::

    int xvel_bc_type(int marker)
      { return (marker != 2) ? BC_ESSENTIAL : BC_NONE; }

Returning BC_NONE for some part of the boundary assigns degrees of freedom but turns
off all surface integral processing on that part of the boundary, which is what we need
in this case.

Next we rewrite the weak formulation so that it fits into the block form :eq:`weaksystem`:

.. math::
    :nowrap:

    \begin{eqnarray*}   a_{11}(u_1, v_1) &=& \int_\Omega \frac{u_1 v_1}{\tau} \dx +                        \int_\Omega \frac{1}{\rm Re} \nabla u_1 \cdot \nabla v_1 \dx +                        \int_\Omega (\bfu^{n-1} \cdot \nabla) u_1 v_1 \dx, \\   a_{22}(u_2, v_2) &=& \int_\Omega \frac{u_2 v_2}{\tau} \dx +                        \int_\Omega \frac{1}{\rm Re} \nabla u_2 \cdot \nabla v_2 \dx +                        \int_\Omega (\bfu^{n-1} \cdot \nabla) u_2 v_2 \dx, \end{eqnarray*}


.. math::
    :nowrap:

    \begin{eqnarray*}   a_{13}(p, v_1) &=& -\int_\Omega p \dd{v_1}{x} \dx, \\   a_{23}(p, v_2) &=& -\int_\Omega p \dd{v_2}{y} \dx, \\   a_{31}(u_1, q) &=&  \int_\Omega \dd{u_1}{x} q \dx, \\   a_{32}(u_2, q) &=&  \int_\Omega \dd{u_2}{y} q \dx, \\   l_1(v_1) &=& \int_\Omega \frac{u^{n-1}_1 v_1}{\tau}, \\   l_2(v_2) &=& \int_\Omega \frac{u^{n-1}_2 v_2}{\tau}. \end{eqnarray*}

Notice first that the forms $a_{11}$ and $a_{22}$ are identical, i.e., $a_{11}(u,v) = a_{22}(u,v)$.
Further, the first two terms of $a_{11}$ and $a_{22}$ are symmetric. We will also exploit the
antisymmetry $a_{13}(u,v) = -a_{31}(u,v)$ and $a_{23}(u,v) = -a_{32}(u,v)$ in the following.

The implementation of the symmetric terms in $a_{11}$ and $a_{22}$ is straightforward. The form
\verb"bilinear_form_sym_0_0_1_1" (the same form is used for both $a_{11}$ and $a_{22}$)
simply contains the command
::

    return int_grad_u_grad_v(fu, fv, ru, rv) / Re +
           int_u_v(fu, fv, ru, rv) / tau;

As for the convection term, we need access to the solution on the previous time level, $\bfu^{n-1}$.
This is accomplished by defining two instances of the class Solution at the global level:
::

    // velocities from the previous time step
    Solution xprev, yprev;

In \verb"bilinear_form_unsym_0_0_1_1", which completes the forms $a_{11}$ and $a_{22}$, we can use
the predefined integral \verb"int_w_nabla_u_v" (see the
file {\tt src/integrals\_h1.h})
and plug in {\tt xprev} and {\tt yprev} for the velocity:
::

    return int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv);

The rest of the forms are easy and will not be discussed here. However, there is one more important
thing you need to do if you use external functions (such as xprev and yprev in the
weak forms. Hermes needs to be told about all such functions and where they are used in the weak
formulation, so that they can be initialized properly and also incorporated in the multi-mesh assembling,
if necessary.
Apart from the symmetry flag and the integration area,
add_biform() takes one more optional argument, the number of external functions used by the form,
followed by that many pointers to the external functions. The complete WeakForm initialization
looks like this:
::

    // set up weak formulation
    WeakForm wf(3);
    wf.add_biform(0, 0, callback(bilinear_form_sym_0_0_1_1), SYM);
    wf.add_biform(0, 0, callback(bilinear_form_unsym_0_0_1_1), UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(1, 1, callback(bilinear_form_sym_0_0_1_1), SYM);
    wf.add_biform(1, 1, callback(bilinear_form_unsym_0_0_1_1), UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(0, 2, callback(bilinear_form_unsym_0_2), ANTISYM);
    wf.add_biform(1, 2, callback(bilinear_form_unsym_1_2), ANTISYM);
    wf.add_liform(0, callback(linear_form), ANY, 1, &xprev);
    wf.add_liform(1, callback(linear_form), ANY, 1, &yprev);

Notice also the use of the ANTISYM flag for the forms $a_{13}$ and $a_{23}$, which
saves us a little assembling time and the need to define $a_{31}$ and $a_{32}$.

Before entering the main iteration loop, we need to initialize the previous solutions
xprev and yprev with the initial condition 
:eq:`ns:initial`. Besides holding the finite element solution, the Solution class
can be forced to return zero, to return a constant, or to return an arbitrary function
using the methods set_zero(), set_const() and set_exact(), respectively
Here we simply call set_zero() and supply the function domain, i.e., the mesh:
::

    // initial BC: xprev and yprev are zero
    xprev.set_zero(&mesh);
    yprev.set_zero(&mesh);

We are now ready to start the iterative process. In each iteration, we assemble the
stiffness matrix and solve for the unknown velocity xsln, ysln and
pressure psln on the current time level:
::

    // assemble and solve
    Solution xsln, ysln, psln;
    sys.assemble();
    sys.solve(3, &xsln, &ysln, &psln);

At the end of each iteration, the current solution must be remembered as the future
previous solution. This is done by assigning xsln and ysln to xprev
and yprev:
::

    xprev = xsln;
    yprev = ysln;

The assignment operator is overloaded for Solution and in fact is equal to calling
Solution::assign(), which is an efficient way of handing over solution data from
one Solution to another.
The velocity is visualized in each iteration using the VectorView class, as shown
in the following figure:

.. image:: img/velocity.jpg
   :align: center
   :width: 600
   :height: 260
   :alt: Velocity solution visualized with the class VectorView.

Waveguide
---------

To be added soon.

Crack
-----

To be added soon.

Nernst-Planck
-------------

To be added soon.

Gross-Pitaevski
---------------

To be added soon.

