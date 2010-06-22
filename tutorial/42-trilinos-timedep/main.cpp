#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  The purpose of this example is to show how to use Trilinos
//  for time-dependent PDE problem.
//  NOX solver is used, either using Newton's method or JFNK and
//  with or without preconditioning,
//
//  PDE: Heat transfer.
//
//  Domain: Unit square.
//
//  BC: Dirichlet at the bottom, Newton elsewhere.
//

const int INIT_REF_NUM = 4;       // Number of initial uniform mesh refinements.
const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const double ALPHA = 10.0;
const double LAMBDA = 1e5;
const double HEATCAP = 1e6;
const double RHO = 3000.0;
const double TEMP_EXT = 20.0;
const double TEMP_INIT = 10.0;

const double TAU = 50.0;          // Time step.        

const bool JFNK = true;
const bool PRECOND = true;

const int BDY_BOTTOM = 1;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == BDY_BOTTOM) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Essential (Dirichlet) boundary conditions values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return TEMP_INIT;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  info("Number of DOF: %d", space.get_num_dofs());

  // Solutions.
  Solution tprev, titer, tsln;

  // Define constant initial condition. 
  tprev.set_const(&mesh, 20.0);
  titer.set_const(&mesh, 20.0);

  // Initialize the weak formulation.
  WeakForm wf(1, JFNK ? true : false);
  wf.add_jacform(callback(jacobian));
  wf.add_jacform_surf(callback(jacobian_surf));
  wf.add_resform(callback(residual), H2D_ANY, &tprev);
  wf.add_resform_surf(callback(residual_surf));

  // Initialize the finite element problem.
  H1Shapeset shapeset;
  FeProblem fep(&wf);
  fep.set_spaces(&space);
  PrecalcShapeset pss(&shapeset);
  //fep.set_pss(1, &pss);

  // Project the function "titer" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solution on the FE mesh.");

  LinSystem ls(&wf, &space);
  ls.project_global(&tprev, &tprev);

  // Get the coefficient vector.
  scalar *vec = ls.get_solution_vector();

  // Measure the projection time.
  double proj_time = cpu_time.tick().last();

  // Initialize NOX solver.
  NoxSolver solver(&fep);

  // Select preconditioner.
  MlPrecond pc("sa");
  if (PRECOND)
  {
    if (JFNK) solver.set_precond(&pc);
    else solver.set_precond("ML");
  }

  // Initialize the view.
  ScalarView Tview("Temperature", 0, 0, 450, 600);
  Tview.set_min_max_range(10,20);

  // Time stepping loop:
  double total_time = 0.0;
  cpu_time.tick_reset();
  for (int ts = 1; total_time <= 2000.0; ts++)
  {
    info("---- Time step %d, t = %g s", ts, total_time += TAU);

    info("Assembling by FeProblem, solving by NOX.");
    solver.set_init_sln(vec);
    bool solved = solver.solve();
    if (solved)
    {
      vec = solver.get_solution_vector();
      tsln.set_fe_solution(&space, &pss, vec);
      Tview.show(&tsln);
      tprev = tsln;
    }
    else
      error("NOX failed.");

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      solver.get_num_iters(), solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      solver.get_num_lin_iters(), solver.get_achieved_tol());
  }

  info("Total running time: %g", cpu_time.accumulated());

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
