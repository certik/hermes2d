=======================================
Tutorial Part II (Automatic Adaptivity)
=======================================

In the computations that we carried out so far, we have not paid any attention
to the accuracy of the results. In general, a computation on a fixed mesh is
not likely to be very accurate. There is a need for *adaptive mesh refinement
(AMR)* algorithms that improve the quality of the approximation by refining
mesh elements where the approximation is bad.

Adaptive h-FEM and hp-FEM
-------------------------

In traditional low-order FEM, refining an element is not algorithmically complicated,
and so the most difficult part is to find out what elements should be
refined. To do this, people employ various techniques ranging from rigorous
guaranteed a-posteriori error estimates to heuristic criteria such as residual
error indicators, error indicators based on steep gradients, etc. Unfortunately,
none of these approaches is suitable for Hermes: The rigorous guaranteed error
estimates only exist for very simple problems, such as linear elliptic PDEs,
and thus they are far from PDE-independent. Heuristic techniques are not
employed in Hermes for the same reason, and moreover since such criteria
lack a transparent relation to the true approximation error.

Adaptive low-order FEM is known to be notoriously inefficient, and practitioners
are rightfully skeptical of it. The reason is illustrated here:

.. image:: img/lshape/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: Typical convergence curves for adaptive linear FEM, quadratic FEM, and *hp*-FEM.

These convergence curves are typical representative examples, confirmed with
many numerical experiments of independent researchers, and supported with
theory. The horizontal axis shows (in linear scale) the number of degrees of freedom
(= size of the stiffness matrix) that increases during automatic adaptivity. The
vertical one shows the approximation error (in logarithmic scale). Note that in all
three cases, the error drops very fast during a short initial phase of the adaptive
computation. However, with both linear and quadratic FEM, the convergence slows
down dramatically as the adaptivity progresses. Note that the low-order FEM
is doomed to such slow convergence by its poor approximation properties -
an excellent adaptivity algorithm cannot improve it (and a bad
algorithm can make it even worse).

In order to obtain fast, usable adaptivity (the green curve), one
has to resort to adaptive *hp*-FEM. The *hp*-FEM takes advantage of two facts:

* Large high-degree elements approximate smooth parts of solution much better than small linear ones. We created the example 'smooth' to illustrate this fact. Check it out, the results are impressive.
* This holds the other way where the solution is not smooth.

Automatic adaptivity in the *hp*-FEM is substantially different from adaptivity
in low-order FEM, since every element can be refined in many different ways.
The following figure shows several refinement candidates for a fourth-order element.

.. image:: img/refinements.png
   :align: center
   :width: 650
   :height: 300
   :alt: Examples of *hp*-refinements.

Due to the large number of refinement options, classical error estimators (that
provide a constant error estimate per element) cannot be used to guide automatic 
*hp*-adaptivity. For this, one needs to know the *shape* of the
approximation error.

In analogy to the most successful adaptive ODE solvers,
Hermes uses a pair of approximations with different orders of accuracy to obtain
this information: *coarse mesh solution* and 
*fine mesh solution*. The initial coarse mesh is read from the mesh file,
and the initial fine mesh is created through its global refinement both in
$h$ and $p$.
The fine mesh solution is the approximation of interest both during the adaptive
process and at the end of computation. The coarse mesh
solution represents its low-order part.

Both these solutions are evolved during the adaptive process
in a PDE-independent manner, based on the discrepancies between global and local
orthogonal projections. (Sometimes we replace the global orthogonal projection with
the solve on the coarse mesh, the difference is negligible.)

The obvious disadvantage of this approach to adaptivity is its higher computational cost,
especially in 3D. We are aware of this fact and would not mind at all replacing it with
some cheaper technique (as long as it also is PDE-independent, works for elements of high 
orders, and can be successfully used to guide *hp*-adaptivity).

Understanding Convergence Rates
-------------------------------

Hermes provides convergence graphs for every adaptive computation. Therefore,
let us spend a short moment explaining their meaning.
The classical notion of $O(h^p)$ convergence rate is related to sequences of 
uniform meshes with a gradually decreasing diameter $h$. In $d$ spatial dimensions, 
the diameter $h$ of a uniform mesh is related to the number of degrees of freedom $N$
through the relation 

.. math::

    h = O(N^{-p/d}).

Therefore a slope of $-p/d$ on the log-log scale means that $err \approx O(N^{-p/d})$
or $err \approx O(h^p)$. When local refinements are enabled, the meaning of $O(h^p)$
convergence rate loses its meaning, and one should switch to convergence in terms of 
the number of degrees of freedom (DOF) or CPU time - Hermes provides both. 

Algebraic convergence of adaptive :math:`h`-FEM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using elements of degree $p$, the convergence rate of adaptive $h$-FEM will not exceed the 
one predicted for uniformly refined meshes (this can be explained using 
mathematical analysis). Nevertheless, the convergence may be faster due to a different 
constant in front of the $h^p$ term. This is illustrated in the following two figures,
both of which are related to a 2D problem with known exact solution. The first pair of 
graphs corresponds to adaptive $h$-FEM with linear elements. The slope on the log-log
graph is -1/2 which means first-order convergence, as predicted by theory. 

.. image:: img/conv-intro/layer_h1.png
   :align: center
   :width: 600
   :height: 450
   :alt: Convergence graph.

The next pair of convergence graphs corresponds to adaptive $h$-FEM with quadratic elements. 
The slope on the log-log graph is -1, which means that the convergence is quadratic as 
predicted by theory.

.. image:: img/conv-intro/layer_h2.png
   :align: center
   :width: 600
   :height: 450
   :alt: Convergence graph.

Note that one always should look at the end of the convergence curve, not at the 
beginning. The automatic adaptivity in Hermes is guided with the so-called 
*reference solution*, which is an approximation on a globally-refined mesh.
In early stages of adaptivity, the reference solution and in turn also the error 
estimate usually are not sufficiently accurate to deliver the expected convergence 
rates. 

Exponential convergence of adaptive $hp$-FEM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is predicted by theory that adaptive $hp$-FEM should attain 
exponential convergence rate. This means that the slope of the
convergence graph is steadily increasing, as shown in the 
following figure.

.. image:: img/conv-intro/aniso-hp.png
   :align: center
   :width: 600
   :height: 450
   :alt: Convergence graph.

While this often is the case with adaptive $hp$-FEM, there are 
problems whose difficulty is such that the convergence is not 
exponential. Or at least not during a long pre-asymptotic 
stage of adaptivity. This may happen, for example, when the solution 
contains an extremely strong singularity. Then basically all error 
is concentrated there, and all adaptive methods will do the same, 
which is to throw into the singularity as many small low-order 
elements as possible. Then the convergence of adaptive $h$-FEM 
and $hp$-FEM may be very similar (usually quite poor).


Estimated vs. exact convergence rates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whenever exact solution is available, Hermes provides both 
estimated error (via the reference solution) as well as the 
exact error. Thus the user can see the quality of the 
error estimate. Note that the estimated error usually is 
slightly less than the exact one, but during adaptivity 
they quickly converge together and become virtually identical. 
This is shown in the figure below.

.. image:: img/conv-intro/layer-hp.png
   :align: center
   :width: 600
   :height: 450
   :alt: Convergence graph.

Electrostatic Micromotor Problem (10)
-------------------------------------

**Git reference:** Tutorial example `10-adapt <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/10-adapt>`_. 

Let us demonstrate the use of adaptive h-FEM and hp-FEM on a linear elliptic problem
concerned with the calculation of
the electrostatic potential in the vicinity of the electrodes of an electrostatic
micromotor. This is a MEMS device free of any coils, and thus resistive to
strong electromagnetic waves (as opposed to classical electromotors).
The following figure shows one half of the domain $\Omega$
(dimensions need to be scaled with $10^{-5}$ and are in meters):

.. image:: img/micromotor.png
   :align: center
   :width: 550
   :height: 400
   :alt: Computational domain for the micromotor problem.

The subdomain $\Omega_2$ represents the moving part of the domain and the area bounded by $\Gamma_2$
represents the electrodes that are fixed. The distribution of the electrostatic potential $\varphi$ is governed by the equation

.. math::

    -\nabla\cdot\left(\epsilon_r\nabla\varphi\right) = 0,

equipped with the Dirichlet boundary conditions

.. math::

    \varphi = 0 V \ \ \ \ \ \mbox{on}\ \Gamma_1,


.. math::

    \varphi = 50 V \ \ \ \ \mbox{on}\ \Gamma_2.

The relative permittivity $\epsilon_r$ is piecewise-constant, $\epsilon_r = 1$ in $\Omega_1$ and
$\epsilon_r = 10$ in $\Omega_2$. The weak formulation reads

.. math::

    \int_\Omega \epsilon_r \nabla u \cdot \nabla v \dx = 0.

The varying parameter $\epsilon_r$ is handled by defining two bilinear forms in the code, one for
$\Omega_1$ and the other for $\Omega_2$. These two areas are delimited by element markers 1 and 2 in
the mesh, and the two forms are assigned to the corresponding markers during the registration of
the forms:

::

    WeakForm wf(1);
    wf.add_biform(0, 0, callback(biform1), H2D_SYM, 1);
    wf.add_biform(0, 0, callback(biform2), H2D_SYM, 2);

The principal part of the example is the main adaptivity loop. However, before the loop is entered, a selector should be created:
::

    H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);

The selector is used by the class H1Adapt to determine how an element should be refined. For that purpose, the selector does following steps:

#. It generates candidates.  A candidate is a proposed refinement.
#. It estimates their local error by projecting the reference solution onto the refinement candidates.
#. It calculates a number of degree of freedom (DOF).
#. It calculates a score of every candidate and sorts candidates according to the score.
#. It selects a candidate with the highest score. If the next candidate has almost the same score, it skips both of them since it is not able to decide which is the best. This helps to keep a mesh symmetric if anisotropic refinements are used on a symmetric solution.

By default, the score is

.. math::

    s = \frac{\log_{10} e_0 - \log_{10} e}{(d_0 - d)^\xi},

where $e$ and $d$ are an estimated error and an estimated number of DOF of a candidate respectively, $e_0$ and $d_0$ are an estimated error and an estimated number of DOF of the examined element respectively, and $\xi$ is a convergence exponent.

The first parameter ``CAND_LIST`` specifies which candidates are generated. In a case of quadrilaterals, all possible values and considered candidates are summarized in the following table:

.. image:: img/cand_list.quads.*
   :align: center
   :alt: Candidates generated for a given candidate list.

The second parameter ``CONV_EXP`` is a convergence exponent used to evaluate a score.

The third parameter specifies the the maximum considered order used in the resulting refinement. In this case, a constant ``H2DRS_DEFAULT_ORDER`` is used. The constant is defined by Hermes2D library and it corresponds to the maximum order supported by the selector. In this case, this is 9.

Furthermore, the selector allows you to weight errors though a method set_error_weights(). Error weights are applied before the error of a candidate is passed to the calculation of the score. Through this method it is possible to set a preference for a given type of a candidate, i.e., H-candidate, P-candidate, and ANISO-candidate. The error weights can be set anytime and setting error weights to appropriate values can lead to a lower number of DOF. However, the best values of weights depend on a solved problem.

In this particular case, a default error weights are used. The default weights prefer the P-candidate and they are defined as:

- H-candidate weight: $2.0$ (see a constant ``H2DRS_DEFAULT_ERR_WEIGHT_H``)
- P-candidate weight: $1.0$ (see a constant ``H2DRS_DEFAULT_ERR_WEIGHT_P``)
- ANISO-candidate weight: $\sqrt{2.0}$ (see a constant ``H2DRS_DEFAULT_ERR_WEIGHT_ANISO``)

Since these weights are default, it is not necessary to express them explicitly. Nevertheless, if expressed, a particular line of the code would be:
::

    selector.set_error_weights(2.0, 1.0, sqrt(2.0));

Besides the error weights, the selector allows you to modify a default behaviour through the method set_option(). The behavior can be modified anytime. Currently, the method accepts following options:

- ``H2D_PREFER_SYMMETRIC_MESH``: Prefer symmetric mesh when selection of the best candidate is done. If set and if two or more candidates has the same score, they are skipped. This option is set by default.
- ``H2D_APPLY_CONV_EXP_DOF``: Use $d^c - d_0^c$, where $c$ is the convergence exponent, instead of $(d - d_0)^c$ to evaluate the score. This options is *not* set by default.

In this case, default settings are used. If expressed explicitly, the code would be:
::

    selector.set_option(H2D_PREFER_SYMMETRIC_MESH, true);
    selector.set_option(H2D_APPLY_CONV_EXP_DOF, false);

After the selector has been created, an adaptivity loop can be entered. The adaptivity loop is an ordinary while-loop or a for-loop. In each iteration of the adaptivity loop, the coarse problem is solved first:
::

    LinSystem ls(&wf, &solver);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln_coarse);

Next, the reference solution is computed on a globally refined copy of the mesh. For this purpose
a temporary space with (optionally) increased element orders is defined and an extra linear system is assembled and solved.
However, in the most of the cases, this can be done using the class RefSystem. All it needs is a pointer
to an instance of the class LinSystem:

::

    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &sln_fine);

The constructor of the class RefSystem allows a user
to choose a different polynomial degree increment (default value 1)
and another element refinement (default value 1) - see the file 
`src/refsystem.h <http://git.hpfem.org/hermes2d.git/blob/HEAD:/src/refsystem.h>`_:

::

    RefSystem(LinSystem* base, int order_increase = 1, int refinement = 1);

In particular, one may want to use order_increase = 0 for h-adaptivity, and 
order_increase = 2 or 3 at the very beginning of computation when the reference 
mesh is still very coarse and thus the reference solution does not give a meaningful 
error estimate. 
 
In the third and last step of each iteration, we refine our mesh and polynomial degrees stored
in our space using a class H1Adapt. This class offers two services:

* it estimates the overall error of the coarse solution in $H^1$ norm,
* it selects elements with the highest error and uses a user-supplied selector to find a refinement.

In a general case, the class H1Adapt is initialized with a list of spaces through the class Tuple. Since we use just a single component in this case, a pointer to an instance of the class H1Space is the only parameter. Then, the user sets the coarse solution and the fine solution and evaluates the error. By default, the error is calculated as

.. math::

    e = \frac{|| u - u_{ref} ||_{H^1}}{|| u_{ref} ||_{H^1}}.

In the code this looks as follows:

::

    H1Adapt hp(&space);
    hp.set_solutions(&sln_coarse, &sln_fine);
    double err_est = hp.calc_error() * 100;

Finally, if ``err_est`` is still above the threshold ``ERR_STOP``, we perform one
adaptivity step:

::

    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      ndof = assign_dofs(&space);
      if (ndof >= NDOF_STOP) done = true;
    }

The constants ``THRESHOLD``, ``STRATEGY`` and ``MESH_REGULARITY``
have the following meaning:

The constant ``STRATEGY`` indicates which adaptive strategy is used. In all cases, the strategy is applied to elements in an order defined through the error. If the user request to process an element outside this order, the element is processed regardless the strategy. Currently, Hermes2D supportes following strategies:

* ``STRATEGY == 0``: Refine elements until sqrt(``THRESHOLD``) times total error is processed. If more elements have similar error refine all to keep the mesh symmetric.
* ``STRATEGY == 1``: Refine all elements whose error is bigger than ``THRESHOLD`` times the error of the first processed element, i.e., the maximum error of an element.
* ``STRATEGY == 2``: Refine all elements whose error is bigger than ``THRESHOLD``.

The constant ``MESH_REGULARITY``
specifies maximum allowed level of hanging nodes: -1 means arbitrary-level
hanging nodes (default), and 1, 2, 3, ... means 1-irregular mesh,
2-irregular mesh, etc. Hermes does not support adaptivity on regular meshes
because of its extremely poor performance.

It is a good idea to spend some time playing with these parameters to
get a feeling for adaptive *hp*-FEM. Also look at other adaptivity examples in
the examples/ directory: layer, lshape deal with elliptic problems and have
known exact solutions. So do examples screen, bessel for time-harmonic
Maxwell's equations. These examples allow you to compare the error estimates
computed by Hermes with the true error. Examples crack, singpert show
how to handle cracks and singularly perturbed problems, respectively. There
are also more advanced examples illustrating automatic adaptivity for nonlinear
problems solved via the Newton's method, adaptive multimesh *hp*-FEM,
adaptivity for time-dependent problems on dynamical meshes, etc.

But let's return to the micromotor example for a moment again: The computation
starts with a very coarse mesh consisting of a few quadrilaterals, some
of which are moreover very ill-shaped. Thanks to the anisotropic refinement
capabilities of the selector, the mesh quickly adapts to the solution
and elements of reasonable shape are created near singularities, which occur
at the corners of the electrode. Initially, all elements of the mesh
are of a low degree, but as the *hp*-adaptive process progresses, the elements
receive different polynomial degrees, depending on the local smoothness of the
solution.

The gradient was visualized using the class VectorView. We have
seen this in the previous section. We plug in the same solution for both vector
components, but specify that its derivatives should be used:

::

    gview.show(&sln, &sln, H2D_EPS_NORMAL, H2D_FN_DX_0, H2D_FN_DY_0);

.. image:: img/motor-sln.png
   :align: left
   :width: 300
   :height: 300
   :alt: Solution - electrostatic potential $\varphi$ (zoomed).

.. image:: img/motor-grad.png
   :align: right
   :width: 300
   :height: 300
   :alt: Gradient of the solution $E = -\nabla\varphi$ and its magnitude (zoomed).

.. raw:: html

   <hr style="clear: both; visibility: hidden;">

.. image:: img/motor-orders.png
   :align: center
   :width: 300
   :height: 300
   :alt: Polynomial orders of elements near singularities (zoomed).

Convergence graphs of adaptive h-FEM with linear elements, h-FEM with quadratic elements
and hp-FEM are shown below.

.. image:: img/example-10/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph for tutorial example 10-adapt.

The following graph shows convergence in terms of CPU time. 

.. image:: img/example-10/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph for tutorial example 10-adapt.

The Multimesh hp-FEM
--------------------

In multiphysics PDE systems (or just PDE systems) it can happen that one
physical field (solution component) has a singularity or a boundary layer 
where other fields are smooth. If one approximates all fields on the 
same mesh, then the necessity to refine the mesh at the singularity
or boundary layer implies new degrees of freedom for the smooth fields 
as well. This can be very wasteful indeed, as we will see in the next
example that deals with a simplified Fitzhugh-Nagumo system. But let us 
first explain briefly the main idea of the multimesh discretization 
method that we developed to circumvent this problem.

Hermes makes it possible to approximate them 
on individual meshes. These meshes are not completely independent
of each other -- they have a common coarse mesh that we call *master mesh*.
The master mesh is there for algorithmic purposes only, it may not 
even be used for discretization purposes: Every mesh in the system 
is obtained from it via an arbitrary sequence of elementary refinements.
This is illustrated in the following figure, where (A) is the master mesh,
(B) - (D) three different meshes (say, for a coupled problem with three
equations), and (E) is the virtual *union mesh* that is used for assembling.

.. image:: img/multimesh/multimesh.png
   :align: center
   :width: 750
   :alt: Multimesh

The union mesh is not constructed physically in the computer memory -- 
merely it serves as a hint to correctly transform integration points
while integrating over sub-elements of the elements of the existing meshes. 
The following figure shows the integration over an element $Q_k$ of the 
virtual union mesh, and what are the appropriate subelements of the 
existing elements where this integration is performed:

.. image:: img/multimesh/multimesh2.png
   :align: center
   :width: 600
   :alt: Multimesh

As a result, the multimesh discretization of the PDE system is *monolithic*
in the sense that *no physics is lost* -- all integrals in the 
discrete weak formulations are evaluated exactly up to the error in the 
numerical quadrature. In particular, we do not perform operator splitting 
or commit errors while transferring solution data between different meshes.
The multimesh assembling in Hermes works with all meshes at the same time, 
there is no such thing as interpolating or projecting functions between 
different meshes. More details about this method can be found in the 
corresponding `research article <http://science.atmoshome.net/science?_ob=MImg&_imagekey=B6TYH-4X1J73B-V-8Y&_cdi=5619&_user=10&_pii=S0377042709005731&_orig=browse&_coverDate=08%2F18%2F2009&_sk=999999999&view=c&wchp=dGLbVzz-zSkWz&md5=6552d3390232dcffc9ca97e9bb626fb0&ie=/sdarticle.pdf>`_. 

Adaptivity in the Multimesh hp-FEM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In principle, the adaptivity procedure for single PDE could be extended 
directly to systems of PDEs. In other words, two spaces can be passed into a constructor of the class H1Adapt,
two coarse and two fine solutions can be passed into set_solutions(),
and finally, calc_error() and adapt() can be called as before. In this way, error estimates in
$H^1$ norm are calculated for elements in both spaces independently and the
worst ones are refined. However, this approach is not optimal if the PDEs are
coupled, since an error caused in one solution component influences the errors
in other components and vice versa.

Recall that in elliptic problems the bilinear form $a(u,v)$ defines the energetic inner product,

.. math::

    (u,v)_e = a(u,v).

The norm induced by this product,

.. math::

    ||u||_e = \sqrt{(u,u)_e},

is called the *energy norm*. When measuring the error in the energy norm
of the entire system, one can reduce the above-mentioned difficulties dramatically.
When calculating the error on an element, the energy norm accounts
also for the error caused by other solution components. 

It is also worth mentioning that the adaptivity algorithm does not make distinctions 
between various meshes. The elements of *all meshes in the system* are put into one
single array, sorted according to their estimated errors, and then the ones with the 
largest error are refined. In other words, it may happen that all elements marked for refinement 
will belong just to one mesh.

If norms of components are substantially different, it is more beneficial to use a relative error of an element rather than an absolute error. The relative error of an element is an absolute error divided by a norm of a component. This behavior can be requested while calling the method calc_error():

::

    hp.calc_error(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL)

The input parameter of the method calc_error() is a combination that is a pair: one member of the pair has to be a constant ```H2D_TOTAL_ERROR_*```, the other member has to be a constant ```H2D_ELEMENT_ERROR_*```. If not specified, the default pair is ```H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_ABS```. Currently available contants are:

- ```H2D_TOTAL_ERROR_REL```: A returned total error will be an absolute error divided by a total norm.
- ```H2D_TOTAL_ERROR_ABS```: A returned total error will be an absolute error.
- ```H2D_TOTAL_ERROR_REL```: An element error which is used to select elements for refinement will be an absolute error divided by a norm of a corresponding component.
- ```H2D_TOTAL_ERROR_ABS```: An element error which is used to select elements for refinement will be an absolute error.


Simplified Fitzhugh-Nagumo System (11)
--------------------------------------

**Git reference:** Tutorial example `11-adapt-system <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/11-adapt-system>`_. 

We consider a simplified version of the Fitzhugh-Nagumo equation.
This equation is a~prominent example of activator-inhibitor systems in two-component reaction-diffusion 
equations, It describes a prototype of an excitable system (e.g., a neuron) and its stationary form 
is

.. math::

    -d^2_u \Delta u - f(u) + \sigma v = g_1,\\
    -d^2_v \Delta v - u + v = g_2.

Here the unknowns $u, v$ are the voltage and $v$-gate, respectively, 
The nonlinear function 

.. math::

    f(u) = \lambda u - u^3 - \kappa
 
describes how an action potential travels through a nerve. Obviously this system is nonlinear.
In order to make it simpler for this tutorial, we replace the function $f(u)$ with just $u$:

.. math::

    f(u) = u.

The original nonlinear version of this example is planned for inclusion in benchmarks. 

Our computational domain is the square $(-1,1)^2$ and we consider zero Dirichlet conditions 
for both $u$ and $v$. In order to enable fair convergence comparisons, we will use the following 
functions as the exact solution:

.. math::

    u(x,y) = \cos\left(\frac{\pi}{2}x\right) \cos\left(\frac{\pi}{2}y\right),\\
    v(x,y) = \hat u(x) \hat u(y)

where

.. math::

    \hat u(x) = 1 - \frac{e^{kx} + e^{-kx}}{e^k + e^{-k}}

is the exact solution of the one-dimensional singularly perturbed 
problem 

.. math::

    -u'' + k^2 u - k^2 = 0

in $(-1,1)$, equipped with zero Dirichlet boundary conditions. The functions $u$ 
and $v$ defined above evidently satisfy the given boundary conditions, and 
they also satisfy the equation, since we inserted them into the PDE system 
and calculated the source functions $g_1$ and $g_2$ from there. These functions 
are not extremely pretty, but they are not too bad either:

::

    // functions g_1 and g_2
    double g_1(double x, double y) 
    {
      return (-cos(M_PI*x/2.)*cos(M_PI*y/2.) + SIGMA*(1. - (exp(K*x) + exp(-K*x))/(exp(K) + exp(-K))) 
             * (1. - (exp(K*y) + exp(-K*y))/(exp(K) + exp(-K))) + pow(M_PI,2.)*pow(D_u,2.)*cos(M_PI*x/2.)
             *cos(M_PI*y/2.)/2.);
    }

    double g_2(double x, double y) 
    {
      return ((1. - (exp(K*x) + exp(-K*x))/(exp(K) + exp(-K)))*(1. - (exp(K*y) + exp(-K*y))/(exp(K) + exp(-K))) 
             - pow(D_v,2.)*(-(1 - (exp(K*x) + exp(-K*x))/(exp(K) + exp(-K)))*(pow(K,2.)*exp(K*y) + pow(K,2.)*exp(-K*y))/(exp(K) + exp(-K)) 
             - (1. - (exp(K*y) + exp(-K*y))/(exp(K) + exp(-K)))*(pow(K,2.)*exp(K*x) + pow(K,2.)*exp(-K*x))/(exp(K) + exp(-K))) - 
             cos(M_PI*x/2.)*cos(M_PI*y/2.));

    }

The weak forms can be found in the 
file `forms.cpp <http://git.hpfem.org/hermes2d.git/blob/HEAD:/tutorial/11-adapt-system/forms.cpp>`_ and 
they are registered as follows:

::

    // initialize the weak formulation
    WeakForm wf(2);
    wf.add_biform(0, 0, callback(bilinear_form_0_0));  
    wf.add_biform(0, 1, callback(bilinear_form_0_1));  
    wf.add_biform(1, 0, callback(bilinear_form_1_0));
    wf.add_biform(1, 1, callback(bilinear_form_1_1));
    wf.add_liform(0, linear_form_0, linear_form_0_ord);
    wf.add_liform(1, linear_form_1, linear_form_1_ord);

Beware that despite each of the forms is actually symmetric, one cannot use the H2D_SYM flag as in the 
elasticity equations, since it has a slightly different 
meaning (see example `08-system <http://hpfem.org/hermes2d/doc/src/tutorial-1.html#systems-of-equations-08>`_).

The exact error is for us the maximum of the relative errors of the two solution components:

::

    // calculate error wrt. exact solution
    ExactSolution uexact(&umesh, u_exact);
    ExactSolution vexact(&vmesh, v_exact);
    double u_error = h1_error(&u_sln_coarse, &uexact) * 100;
    double v_error = h1_error(&v_sln_coarse, &vexact) * 100;
    double error = fmax(u_error, v_error);
    info("Exact solution error for u (H1 norm): %g%%", u_error);
    info("Exact solution error for v (H1 norm): %g%%", v_error);
    info("Exact solution error (maximum): %g%%", error);

The next code snippet shows how we define the energy norm for adaptive
multimesh hp-FEM, whose necessity was explained in the previous 
paragraph:

::

    // calculate element error estimates and the total error estimate
    H1Adapt hp(Tuple<Space*>(&uspace, &vspace));
    hp.set_solutions(Tuple<Solution*>(&u_sln_coarse, &v_sln_coarse), Tuple<Solution*>(&u_sln_fine, &v_sln_fine));
    hp.set_biform(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
    hp.set_biform(0, 1, bilinear_form_0_1<scalar, scalar>, bilinear_form_0_1<Ord, Ord>);
    hp.set_biform(1, 0, bilinear_form_1_0<scalar, scalar>, bilinear_form_1_0<Ord, Ord>);
    hp.set_biform(1, 1, bilinear_form_1_1<scalar, scalar>, bilinear_form_1_1<Ord, Ord>);
    double err_est = hp.calc_error(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_ABS) * 100;
    info("Estimate of error wrt. ref. solution (energy norm): %g%%", err_est);

Notice that in this case, the method calc_error() specifies explicitly that an error of an element is not devided by the norm of a component and the method should return a total error which is relative, i.e., divided by the norm. This is a default behavior.

The following two figures show the solutions $u$ and $v$. Notice their 
large qualitative differences: While $u$ is smooth in the entire domain, 
$v$ has a thin boundary layer along the boundary:

.. image:: img/example-11/solution_u.png
   :align: center
   :width: 465
   :height: 400
   :alt: Solution

.. image:: img/example-11/solution_v.png
   :align: center
   :width: 465
   :height: 400
   :alt: Solution

Resulting mesh for $u$ and $v$ obtained using conventional (single-mesh) hp-FEM: 12026 DOF
(6013 for each solution). 

.. image:: img/example-11/mesh_single.png
   :align: center
   :width: 465
   :height: 400
   :alt: Mesh

Resulting mesh for $u$ obtained using the multimesh hp-FEM: 169 DOF

.. image:: img/example-11/mesh_multi_u.png
   :align: center
   :width: 465
   :height: 400
   :alt: Mesh

Resulting mesh for $v$ obtained using the multimesh hp-FEM: 3565 DOF

.. image:: img/example-11/mesh_multi_v.png
   :align: center
   :width: 465
   :height: 400
   :alt: Mesh

DOF convergence graphs:

.. image:: img/example-11/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: img/example-11/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

Adaptivity for General 2nd-Order Linear Equation (12)
-----------------------------------------------------

**Git reference:** Tutorial example `12-adapt-general <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/12-adapt-general>`_. 

This example does not bring anything substantially new and its purpose is solely to 
save you work adding adaptivity to the tutorial example 
`07-general <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/07-general>`_. 
Feel free to adjust the 
`main.cpp <http://git.hpfem.org/hermes2d.git/blob/HEAD:/tutorial/12-adapt-general/main.cpp>`_ 
file for your own applications.

Solution:

.. image:: img/example-12/12-solution.png
   :align: center
   :width: 465
   :height: 400
   :alt: Solution to the general 2nd-order linear equation example.

Final hp-mesh:

.. image:: img/example-12/12-mesh.png
   :align: center
   :width: 450
   :height: 400
   :alt: Final finite element mesh for the general 2nd-order linear equation example.

Convergence graphs of adaptive h-FEM with linear elements, h-FEM with quadratic elements
and hp-FEM.

.. image:: img/example-12/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph for tutorial example 12-adapt-general.

Convergence comparison in terms of CPU time. 

.. image:: img/example-12/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph for tutorial example 12-adapt-general.

