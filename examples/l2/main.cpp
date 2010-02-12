#include "hermes2d.h"
#include "solver_umfpack.h"

// This example shows how to use the L2 finite element space and L2 shapeset.
// As a sample problem, a continuous function x^3 + y^3 is projected onto the
// L2 finite element space in the L2 norm. When zero-order is used, the result
// is a piecewice constant function.
//
// The following parameters can be changed:

const int P_INIT = 1;

// projected function
double F(double x, double y)
{
  return x*x*x + y*y*y;
}

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
  for (int i = 0; i < n; i++)
    result += wt[i] * ((pow(e->x[i], 3) + pow(e->y[i], 3)) * v->val[i]);
  return result;
}

// boundary conditions
int bc_types(int marker)
{
   return BC_NONE;
}

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
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  // initialize the shapeset and the cache
  L2Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create the L2 space
  L2Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);

  // set uniform polynomial degrees
  space.set_uniform_order(P_INIT);

  // enumerate basis functions
  space.assign_dofs();

/*  BaseView bview;
  bview.show(&space);
  bview.wait_for_close();
*/

  Solution sln;

  // matrix solver
  UmfpackSolver umfpack;

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
  view1.wait_for_keypress();

  // wait for keyboard or mouse input
  View::wait("Waiting for all views to be closed.");
  return 0;
}

