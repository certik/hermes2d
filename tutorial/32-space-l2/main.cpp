#include "hermes2d.h"
#include "solver_umfpack.h"
#include "shapeset_common.h"

// This example shows how to use the L2 finite element space and L2 shapeset.
// As a sample problem, a continuous function x^3 + y^3 is projected onto the
// L2 finite element space in the L2 norm. When zero-order is used, the result
// is a piecewice constant function. The class BaseView will show you the basis
// functions.
//
// The following parameters can be changed:

const int INIT_REF_NUM = 1;    // Number of initial uniform mesh refinements.
const int P_INIT = 3;          // Polynomial degree of mesh elements.

// Projected function.
double F(double x, double y)
{
  return - pow(x, 4) * pow(y, 5); 
}

// Bilinear and linear form defining the projection.
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_u_v<Real, Scalar>(n, wt, u, v);
}

// Returns the value \int v dx.
template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (-pow(e->x[i], 4) * pow(e->y[i], 5)) * v->val[i];
  return result;
}

// Boundary condition types.
BCType bc_types(int marker)
{
   return BC_NONE;
}

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform uniform mesh refinements.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize the shapeset.
  L2Shapeset shapeset;

  // Create an L2 space.
  L2Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);

  // Set uniform polynomial degrees.
  space.set_uniform_order(P_INIT);

  // Enumerate basis functions.
  int ndof = assign_dofs(&space);

  // View basis functions.
  BaseView bview;
  bview.show(&space);
  View::wait(H2DV_WAIT_KEYPRESS);

  // Initialize solution.
  Solution sln;

  // Matrix solver.
  UmfpackSolver solver;

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_biform(callback(bilinear_form));
  wf.add_liform(callback(linear_form));

  // Assemble and solve the finite element problem.
  LinSystem sys(&wf, &solver, &space);
  sys.assemble();
  sys.solve(&sln);

  // Visualize the solution.
  ScalarView view1("Solution 1");
  view1.show(&sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

