#include "hermes2d.h"
#include "solver_umfpack.h"
#include "function.h"

//  This test makes sure that example 13-newton-elliptic-1 works correctly.

const int P_INIT = 2;             // Initial polynomial degree
const double NEWTON_TOL = 1e-6;   // Stopping criterion for the Newton's method
const int INIT_GLOB_REF_NUM = 3;  // Number of initial uniform mesh refinements
const int INIT_BDY_REF_NUM = 5;   // Number of initial refinements towards boundary

// Thermal conductivity (temperature-dependent)
// Note: for any u, this function has to be positive
template<typename Real>
Real lam(Real u) { return 1 + pow(u, 4); }

// Derivative of the thermal conductivity with respect to 'u'
template<typename Real>
Real dlam_du(Real u) { return 4*pow(u, 3); }

// Boundary condition type (essential = Dirichlet)
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Heat sources (can be a general function of 'x' and 'y')
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Jacobian matrix
template<typename Real, typename Scalar>
Scalar jac(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (dlam_du(u_prev->val[i]) * u->val[i] * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
                       + lam(u_prev->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
                       
  return result;
}

// Fesidual vector
template<typename Real, typename Scalar>
Scalar res(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (lam(u_prev->val[i]) * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
		       - heat_src(e->x[i], e->y[i]) * v->val[i]);
  return result;
}

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // initial mesh refinements
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1,INIT_BDY_REF_NUM);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_uniform_order(P_INIT);
  space.assign_dofs();

  // previous solution for the Newton's iteration
  Solution u_prev;

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(jac), UNSYM, ANY, 1, &u_prev);
  wf.add_liform(0, callback(res), ANY, 1, &u_prev);

  // initialize the nonlinear system and solver
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  // set zero function as the initial condition
  u_prev.set_zero(&mesh);
  nls.set_ic(&u_prev, &u_prev);

  // Newton's loop
  int it = 1;
  double res_l2_norm;
  Solution sln;
  do
  {
    info("\n---- Newton iter %d ---------------------------------\n", it++);

    // assemble the Jacobian matrix and residual vector, 
    // solve the system
    nls.assemble();
    nls.solve(1, &sln);

    // calculate the l2-norm of residual vector
    res_l2_norm = nls.get_residuum_l2_norm();
    info("Residuum L2 norm: %g\n", res_l2_norm);

    // save the new solution as "previous" for the 
    // next Newton's iteration
    u_prev = sln;

    if (it > 18) break;
  }
  while (res_l2_norm > NEWTON_TOL);

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  if (it <= 18 && res_l2_norm < 1e-10) {  // it should be 18 and res_l2_norm = 4.87539e-12
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}

