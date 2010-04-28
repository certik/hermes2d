#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"

//  This is an elliptic problem with known exact solution. The solution contains
//  a very strong singularity that represents a challenge for most adaptive methods.
//
//  PDE: -div(A(x,y) grad u) = 0
//  where a(x,y) = R in the first and third quadrant
//               = 1 in the second and fourth quadrant
//
//  Exact solution: u(x,y) = cos(M_PI*y/2)    for x < 0
//                  u(x,y) = cos(M_PI*y/2) + pow(x, alpha)   for x > 0   where alpha > 0
//
//  Domain: Square (-1,1)^2.
//
//  BC:  Dirichlet given by exact solution.

const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;       // number of initial mesh refinements (the original mesh is just one element)
const double THRESHOLD = 0.3;     // This is a quantitative parameter of the adapt(...) function and
                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;           // Adaptive strategy:
                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                  //   error is processed. If more elements have similar errors, refine
                                  //   all to keep the mesh symmetric.
                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                  //   than THRESHOLD times maximum element error.
                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                  //   than THRESHOLD.
                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const RefinementSelectors::AllowedCandidates ADAPT_TYPE = RefinementSelectors::H2DRS_CAND_HP;         // Type of automatic adaptivity.
const bool ISO_ONLY = false;      // Isotropic refinement flag (concerns quadrilateral elements only).
                                  // ISO_ONLY = false ... anisotropic refinement of quad elements
                                  // is allowed (default),
                                  // ISO_ONLY = true ... only isotropic refinements of quad elements
                                  // are allowed.
const int MESH_REGULARITY = -1;   // Maximum allowed level of hanging nodes:
                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                  // Note that regular meshes are not supported, this is due to
                                  // their notoriously bad performance.
const double ERR_STOP = 1.0;      // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;     // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

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

static double fndd(double x, double y, double& dx, double& dy)
{
  double theta = atan2(y,x);
  if (theta < 0) theta = theta + 2*M_PI;
  double r = sqrt(x*x + y*y);
  // x-derivative
  if (theta <= M_PI/2.) {
    dx = TAU*x*pow(r, (2.*(-1 + TAU/2.))) *
    cos((M_PI/2. - SIGMA)*TAU) *
    cos(TAU*(-M_PI/2. + RHO + theta)) +
    (TAU*y*pow(r, TAU)*cos((M_PI/2. - SIGMA)*TAU) *
    sin(TAU*(-M_PI/2. + RHO + theta))/(r*r));
  }
  else {
    if (theta <= M_PI) {
      dx = TAU*x * pow(r, (2.*(-1 + TAU/2.))) * cos(RHO*TAU) *
      cos(TAU*(-M_PI + SIGMA + theta)) +
      (TAU*y * pow(r, TAU) * cos(RHO*TAU) *
      sin(TAU*(-M_PI + SIGMA + theta))/(r*r));
    }
    else {
      if (theta <= 3.*M_PI/2.) {
        dx = TAU*x * pow(r, (2.*(-1 + TAU/2.))) * cos(SIGMA*TAU) *
        cos(TAU*(-M_PI - RHO + theta)) +
	(TAU*y * pow(r, TAU) * cos(SIGMA*TAU) *
	 sin(TAU*(-M_PI - RHO + theta))/(r*r));
      }
      else {
        dx = TAU*x* pow(r, (2*(-1 + TAU/2.))) *
        cos((M_PI/2. - RHO)*TAU) *
	cos(TAU*(-3.*M_PI/2. - SIGMA + theta)) +
	(TAU*y*pow(r, TAU) * cos((M_PI/2. - RHO)*TAU) *
	sin(TAU*(-3.*M_PI/2. - SIGMA + theta))/(r*r));
      }
    }
  }
  // y-derivative
  if (theta <= M_PI/2.) {
    dy = TAU*y * pow(r, (2*(-1 + TAU/2.))) *
      cos((M_PI/2. - SIGMA)*TAU) *
      cos(TAU*(-M_PI/2. + RHO + theta)) -
      (TAU * pow(r, TAU) * cos((M_PI/2. - SIGMA)*TAU) *
       sin(TAU*(-M_PI/2. + RHO + theta))*x/(r*r));
  }
  else {
    if (theta <= M_PI) {
      dy = TAU*y* pow(r, (2*(-1 + TAU/2.))) * cos(RHO*TAU) *
	cos(TAU*(-M_PI + SIGMA + theta)) -
        (TAU * pow(r, TAU) * cos(RHO*TAU) *
	 sin(TAU*(-M_PI + SIGMA + theta))*x/(r*r));
    }
    else {
      if (theta <= 3.*M_PI/2.) {
        dy = TAU*y * pow(r, (2*(-1 + TAU/2.))) * cos(SIGMA*TAU) *
	  cos(TAU*(-M_PI - RHO + theta)) -
	  (TAU * pow(r, TAU) * cos(SIGMA*TAU) *
	   sin(TAU*(-M_PI - RHO + theta))*x/(r*r));
      }
      else {
        dy = TAU*y * pow(r, (2*(-1 + TAU/2.))) *
        cos((M_PI/2. - RHO)*TAU) *
	  cos(TAU*(-3.*M_PI/2. - SIGMA + theta)) -
	  (TAU * pow(r, TAU) * cos((M_PI/2. - RHO)*TAU) *
	   sin(TAU*((-3.*M_PI)/2. - SIGMA + theta))*x/(r*r));
      }
    }
  }
  return fn(x,y);
}

// Boundary condition types
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Boundary condition values
scalar bc_values(int marker, double x, double y)
{
  return fn(x, y);
}

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

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square_quad.mesh", &mesh);

  // initial mesh refinement
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create finite element space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // enumerate degrees of freedom
  int ndof = assign_dofs(&space);

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form_I_III), H2D_SYM, 0);
  wf.add_biform(0, 0, callback(bilinear_form_II_IV), H2D_SYM, 1);

  // visualize solution and mesh
  ScalarView sview("Coarse solution", 0, 100, 798, 700);
  OrderView  oview("Polynomial orders", 800, 100, 798, 700);

  // matrix solver
  UmfpackSolver solver;

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof_est, graph_dof_exact, graph_cpu_est, graph_cpu_exact;

  // prepare selector
  RefinementSelectors::H1NonUniformHP selector(ISO_ONLY, ADAPT_TYPE, 1.0, H2DRS_DEFAULT_ORDER, &shapeset);

  // adaptivity loop
  int it = 1;
  bool done = false;
  TimePeriod cpu_time;
  Solution sln_coarse, sln_fine;
  do
  {
    info("!---- Adaptivity step %d ---------------------------------------------", it); it++;

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the coarse mesh problem
    LinSystem ls(&wf, &solver);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln_coarse);

    // time measurement
    cpu_time.tick();

    // calculate error wrt. exact solution
    ExactSolution exact(&mesh, fndd);
    double error = h1_error(&sln_coarse, &exact) * 100;

    // view the solution
    sview.show(&sln_coarse);
    oview.show(&space);

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the fine mesh problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &sln_fine);

    // calculate error estimate wrt. fine mesh solution
    H1AdaptHP hp(1, &space);
    double err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;

    // time measurement
    cpu_time.tick();

    // report results
    info("Exact solution error: %g%%", error);
    info("Estimate of error: %g%%", err_est);

    // add entries to DOF convergence graphs
    graph_dof_exact.add_values(space.get_num_dofs(), error);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_dof_est.add_values(space.get_num_dofs(), err_est);
    graph_dof_est.save("conv_dof_est.dat");

    // add entries to CPU convergence graphs
    graph_cpu_exact.add_values(cpu_time.accumulated(), error);
    graph_cpu_exact.save("conv_cpu_exact.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est);
    graph_cpu_est.save("conv_cpu_est.dat");

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(THRESHOLD, STRATEGY, &selector, MESH_REGULARITY);
      ndof = assign_dofs(&space);
      if (ndof >= NDOF_STOP) done = true;
    }

    // time measurement
    cpu_time.tick();
  }
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // show the fine solution - this is the final result
  sview.set_title("Final solution");
  sview.show(&sln_fine);

  // wait for all views to be closed
  View::wait();
  return 0;
}
