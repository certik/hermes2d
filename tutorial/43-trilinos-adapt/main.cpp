#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  The purpose of this example is to show how to use Trilinos while adapting mesh
//  Solved by NOX solver, either using Newton's method or JFNK, with or without preconditioning.
//  The underlying problem is the "layer" benchmark with known exact solution.
//
//  PDE: -Laplace u = f.
//
//  Domain: Unit square.
//
//  BC: Nonhomogeneous Dirichlet.
//
//  Known exact solution, see functions fn() and fndd().
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 4;              // Number of initial uniform mesh refinements.
const int P_INIT = 2;                    // Initial polynomial degree of all mesh elements.
const bool JFNK = true;                  // true = jacobian-free method,
                                         // false = Newton.
const bool PRECOND = true;               // Preconditioning by jacobian in case of jfnk,
                                         // default ML proconditioner in case of Newton.
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1.0;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.

// Problem parameters.
double SLOPE = 60;                       // Slope of the layer inside the domain

// Exact solution.
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

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}
// Essential (Dirichlet) boundary conditions.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return fn(x, y);
}

// Right-hand side.
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

// Preconditioner weak form.
template<typename Real, typename Scalar>
Scalar precond_form(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, vi, vj);
}

// Residual weak form.
template<typename Real, typename Scalar>
Scalar residual_form(int n, double *wt, Func<Real> *u[], Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u[0], vj) + int_F_v<Real, Scalar>(n, wt, rhs, vj, e);
}

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  info("Number of DOF: %d", space.get_num_dofs());

  // Initialize the weak formulation.
  WeakForm wf(1, JFNK ? true : false);
  if (PRECOND) wf.add_jacform(callback(precond_form), H2D_SYM);
  wf.add_resform(callback(residual_form));

  // Initialize views.
  ScalarView sview("Coarse mesh solution", 0, 100, 798, 700);
  OrderView  oview("Coarse mesh", 800, 100, 798, 700);

  // DOF convergence graphs.
  SimpleGraph graph_dof_est, graph_dof_exact;

  // Initialize the coarse mesh problem.
  LinSystem ls(&wf, &space);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  double err_exact;
  Solution sln_coarse, sln_fine;
  double *s;
  do
  {
    info("---- Adaptivity step %d:", as);
   
    // Initialize finite element problem.
    H1Shapeset shapeset;
    FeProblem fep(&wf);
    fep.set_spaces(&space);
    PrecalcShapeset pss(&shapeset);
    //fep.set_pss(1, &pss);

    // Initialize NOX solver.
    NoxSolver solver(&fep);

    // Choose preconditioner.
    MlPrecond pc("sa");
    if (PRECOND)
    {
      if (JFNK) solver.set_precond(&pc);
      else solver.set_precond("ML");
    }

    // Solve the matrix problem.
    info("Coarse mesh problem: Assembling by FeProblem, solving by NOX.");
    bool solved = solver.solve();
    if (solved)
    {
      s = solver.get_solution_vector();
      sln_coarse.set_fe_solution(&space, &pss, s);

      info("Coarse Solution info:");
      info(" Number of nonlin iterations: %d (norm of residual: %g)", 
        solver.get_num_iters(), solver.get_residual());
      info(" Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
        solver.get_num_lin_iters(), solver.get_achieved_tol());

      // Time measurement.
      cpu_time.tick();

      // View the solution and mesh.
      sview.show(&sln_coarse);
      oview.show(&space);

      // Skip visualization time.
      cpu_time.tick(H2D_SKIP);
    }

    // FIXME: RefSystem should be used instead of the following:
        // Create uniformly refined reference mesh.
        Mesh rmesh; rmesh.copy(&mesh); 
        rmesh.refine_all_elements();
        // Reference FE space.
        H1Space rspace(&rmesh, bc_types, essential_bc_values, 1);
        int order_increase = 1;
        rspace.copy_orders(&space, order_increase); // increase orders by one
        // Initialize FE problem on reference mesh.
        FeProblem ref_fep(&wf);
        ref_fep.set_spaces(&rspace);
        //ref_fep.set_pss(1, &pss);

    // Initialize NOX solver.
    NoxSolver ref_solver(&ref_fep);
    if (PRECOND)
    {
      if (JFNK) ref_solver.set_precond(&pc);
      else ref_solver.set_precond("ML");
    }

    // Solve the matrix problem using NOX.
    info("Fine mesh problem: Assembling by FeProblem, solving by NOX.");
    solved = ref_solver.solve();
    if (solved)
    {
      s = ref_solver.get_solution_vector();
      sln_fine.set_fe_solution(&rspace, &pss, s);

      info("Reference solution info:");
      info(" Number of nonlin iterations: %d (norm of residual: %g)",
            ref_solver.get_num_iters(), ref_solver.get_residual());
      info(" Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)",
            ref_solver.get_num_lin_iters(), ref_solver.get_achieved_tol());
    }
    else
      error("NOX failed.");

    // Calculate error estimate wrt. fine mesh solution.
    info("Calculating error.");
    H1Adapt hp(&ls);
    hp.set_solutions(&sln_coarse, &sln_fine);
    double err_est = hp.calc_error() * 100;
    ExactSolution exact(&mesh, fndd);
    err_exact = h1_error(&sln_coarse, &exact) * 100;
    info("ndof_coarse: %d, ndof_fine: %d, err_est: %g%%, err_exact: %g%%", 
	 space.get_num_dofs(), rspace.get_num_dofs(), err_est, err_exact);

    // Add entries to DOF convergence graphs.
    graph_dof_exact.add_values(space.get_num_dofs(), err_exact);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_dof_est.add_values(space.get_num_dofs(), err_est);
    graph_dof_est.save("conv_dof_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP) done = true;
    else {
      info("Adapting the coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (space.get_num_dofs() >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (done == false);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
