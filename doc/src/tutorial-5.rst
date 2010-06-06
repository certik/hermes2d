==========================================
Tutorial Part V (Miscellaneous Techniques)
==========================================

This section is a collection of various examples and techniques 
that are worth showing even though they do not fit exactly into 
the previous sections.  

Space H(curl) (30)
------------------

**Git reference:** Tutorial example `30-space-hcurl <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/30-space-hcurl>`_. 

In `example 02-space <http://hpfem.org/hermes2d/doc/src/tutorial-1.html#setting-up-finite-element-space>`_ we first saw how a finite element space over a mesh is created. That was an $H^1$ space suitable for continuous approximations. Another widely used Sobolev space, H(curl), is typically present in Maxwell's problems of electromagnetics. H(curl) approximations are discontinuous, elementwise polynomial vector fields that behave like gradients of $H^1$ functions. (Recall that in electrostatics $E = - \nabla \varphi$.) In particular, H(curl) functions have continuous tangential components along all mesh edges. For the application of the H(curl) space check examples related to Maxwell's equations in the previous sections. Below is a simple code that shows how to set up an H(curl) space and visualize its finite element basis functions:

::

    int P_INIT = 3;

    int main(int argc, char* argv[])
    {
      if (argc < 2) error("Missing mesh file name parameter.");

      // load the mesh
      Mesh mesh;
      H2DReader mloader;
      mloader.load(argv[1], &mesh);

      // uniform mesh refinements
      mesh.refine_all_elements();
      mesh.refine_all_elements();

      // initialize the shapeset and the cache
      HcurlShapeset shapeset;

      // create the Hdiv space
      HcurlSpace space(&mesh, &shapeset);

      // set uniform polynomial degrees
      space.set_uniform_order(P_INIT);

      // enumerate basis functions
      int ndof = assign_dofs(&space);

      // visualise the FE basis
      VectorBaseView bview;
      bview.show(&space);

      // wait for all views to be closed
      View::wait();
      return 0;
    } 

The class VectorBaseView allows the user to browse through 
the finite element basis functions using the left and right 
arrows. A few 
sample basis functions (higher-order bubble functions) are 
shown below. The color shows magnitude of the vector field, 
arrows show its direction.

.. image:: img/example-30/fn0.png
   :align: center
   :width: 300
   :alt: Sample basis function

.. image:: img/example-30/fn1.png
   :align: center
   :width: 300
   :alt: Sample basis function

.. image:: img/example-30/fn2.png
   :align: center
   :width: 300
   :alt: Sample basis function

.. image:: img/example-30/fn3.png
   :align: center
   :width: 300
   :alt: Sample basis function

The space H(curl) is implemented for both quadrilateral and triangular 
elements, and both elements types can be combined in one mesh. 

Space H(div) (31)
-----------------

**Git reference:** Tutorial example `31-space-hdiv <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/31-space-hdiv>`_. 

The space H(div) in 2D is very similar in nature to the space H(curl), except its functions 
behave like (vector-valued) divergences of $H^1$ functions. Finite element basis functions 
in the space H(div) are discontinuous across element interfaces but their normal components 
are continuous. The following code shows how to set up an H(div) space and visualize
its basis functions: 

::

    int P_INIT = 3;

    int main(int argc, char* argv[])
    {
      if (argc < 2) error("Missing mesh file name parameter.");

      // load the mesh
      Mesh mesh;
      H2DReader mloader;
      mloader.load(argv[1], &mesh);

      // uniform mesh refinements
      mesh.refine_all_elements();
      mesh.refine_all_elements();

      // initialize the shapeset and the cache
      HdivShapeset shapeset;

      // create the Hdiv space
      HdivSpace space(&mesh, &shapeset);

      // set uniform polynomial degrees
      space.set_uniform_order(P_INIT);

      // enumerate basis functions
      int ndof = assign_dofs(&space);

      // visualise the FE basis
      VectorBaseView bview;
      bview.show(&space);

      // wait for all views to be closed
      View::wait();
      return 0;
    }

Sample edge functions of polynomial degrees 1, 2, 3, and 4 
corresponding to a boundary edge are shown below:

.. image:: img/example-31/fn0.png
   :align: center
   :width: 300
   :alt: Sample basis function

.. image:: img/example-31/fn1.png
   :align: center
   :width: 300
   :alt: Sample basis function

.. image:: img/example-31/fn2.png
   :align: center
   :width: 300
   :alt: Sample basis function

.. image:: img/example-31/fn3.png
   :align: center
   :width: 300
   :alt: Sample basis function

So far the space H(div) only can be used with quadrilateral elements.

Space L2 (32)
-------------

**Git reference:** Tutorial example `32-space-l2 <http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/31-space-l2>`_. 

We already saw the $L^2$ space in the `Navier-Stokes example <http://hpfem.org/hermes2d/doc/src/tutorial-3.html#navier-stokes-equations>`_ where it was used for pressure to keep the velocity discreetely divergence-free. This example shows how to create an $L^2$ space, visualize 
finite element basis functions, and perform an orthogonal $L^2$-projection of a continuous function onto the FE space. The projected function has the form

::

    // projected function
    double F(double x, double y)
    {
      return x*x*x + y*y*y;
    }

The orthogonal projection is defined via a bilinear form (just an $L^2$ product 
of basis functions) and a linear form ($L^2$ product of basis functions with the 
projected function):

::

    // bilinear and linear form defining the projection
    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_u_v<Real, Scalar>(n, wt, u, v);
    }

    // return the value \int v dx
    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
      result += wt[i] * ((pow(e->x[i], 3) + pow(e->y[i], 3)) * v->val[i]);
      }
      return result;
    }

Here is how to create the space, set a uniform poly degree, enumerate 
degrees of freedom, and show the FE basis:

::

    // create the L2 space
    L2Space space(&mesh, &shapeset);
    space.set_bc_types(bc_types);

    // set uniform polynomial degrees
    space.set_uniform_order(P_INIT);

    // enumerate basis functions
    int ndof = assign_dofs(&space);

    BaseView bview;
    bview.show(&space);
    View::wait(H2DV_WAIT_KEYPRESS);

Next we register the weak forms, assemble and solve 
the matrix problem, and visualize the solution:

::

    // initialize the weak formulation
    WeakForm wf(1);
    wf.add_biform(0, 0, callback(bilinear_form));
    wf.add_liform(0, callback(linear_form));

    // assemble and solve the finite element problem
    LinSystem sys(&wf, &umfpack);
    sys.set_spaces(1, &space);
    sys.set_pss(1, &pss);
    sys.assemble();
    sys.solve(1, &sln);

    // visualize the solution
    ScalarView view1("Solution 1");
    view1.show(&sln);

Sample basis functions:

.. image:: img/example-32/fn0.png
   :align: center
   :width: 400
   :alt: Sample basis function

.. image:: img/example-32/fn1.png
   :align: center
   :width: 400
   :alt: Sample basis function

.. image:: img/example-32/fn2.png
   :align: center
   :width: 400
   :alt: Sample basis function

.. image:: img/example-32/fn3.png
   :align: center
   :width: 400
   :alt: Sample basis function

The projection. Note that this is a discontinuous function:

.. image:: img/example-32/sol.png
   :align: center
   :width: 400
   :alt: Projection

