#include "hermes2d.h"
#include "solver_umfpack.h"
#include "function.h"

//
//  The purpose of this example is to show how to use the Newton's
//  method for a stationary nonlinear PDE problem. Some problem
//  parameters can be changed below.
//
//  PDE: stationary heat transfer with nonlinear thermal conductivity
//  - div[lambda(T)grad T] = 0
//
//  Domain: square
//
//  BC:  T = 100 on the left, top and bottom edges
//       dT/dn = 0 on the right edge
//

/********** Problem parameters ***********/

int P_INIT = 2;            // initial polynomial degree
int PROJ_TYPE = 0;         // 1 for H1 projections, 0 for L2 projections
double NEWTON_TOL = 1e-3;  // convergence criterion for the Newton's method

// thermal conductivity (temperature-dependent)
// for any u, this function has to be  positive in the entire domain!
template<typename Real>
Real lam(Real T) { return 10 + 0.1*pow(T, 2); }
template<typename Real>
Real dlam_dT(Real T) { return 0.1*2*pow(T, 1); }

/********** Definition of boundary conditions ***********/

int bc_types(int marker)
{
 if (marker == 4 || marker == 1 || marker == 3) return BC_ESSENTIAL;
//   if (marker == 4) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

scalar bc_values(int marker, double x, double y)
{
 return 100;
}

/********** Definition of Jacobian matrices and residual vectors ***********/

template<typename Real, typename Scalar>
Scalar residuum(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* titer = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (lam(titer->val[i]) * (titer->dx[i] * v->dx[i] + titer->dy[i] * v->dy[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar jacobian(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* titer = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (dlam_dT(titer->val[i]) * u->val[i] * (titer->dx[i] * v->dx[i] + titer->dy[i] * v->dy[i]) +
                       lam(titer->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  return result;
}

/******************************************************************************/

int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("square.mesh");
  for(int i = 0; i < 5; i++) mesh.refine_all_elements();

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  space.assign_dofs();

  Solution Titer;

  WeakForm wf(1);
  wf.add_biform(0, 0, callback(jacobian), UNSYM, ANY, 1, &Titer);
  wf.add_liform(0, callback(residuum), ANY, 1, &Titer);

  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  char title[100];
  ScalarView view("", 0, 0, 600, 600);
  ScalarView view2("", 700, 0, 600, 600);

  // setting the Dirichlet lift to be the initial condition
  Titer.set_dirichlet_lift(&space, &pss);
  nls.set_ic(&Titer, &Titer, PROJ_TYPE);

  // view initial guess for Newton's method
  /*
  sprintf(title, "Initial guess for the Newton's method");
  view.set_title(title);
  view.show(&Titer);
  view.wait_for_keypress();
  */
  Solution sln;

  int it = 1;
  double res_l2_norm;
  do
  {
    info("\n---- Newton iter %d ---------------------------------\n", it++);

    nls.assemble();
    nls.solve(1, &sln);

    res_l2_norm = nls.get_residuum_l2_norm();
    info("Residuum L2 norm: %g\n", res_l2_norm);
    // want to see Newtons iterations
    sprintf(title, "Newton iteration %d", it-1);
    view.set_title(title);
    view.show(&sln);
    printf("Click into the image window and press any key to proceed.\n");
    view.wait_for_keypress();

    Titer = sln;

  }
  while (res_l2_norm > NEWTON_TOL);

  // wait for keyboard or mouse input
  View::wait("Waiting for keyboard or mouse input.");
  return 0;
}

