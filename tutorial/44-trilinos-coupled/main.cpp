#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

//  The purpose of this example is to show how to use Trilinos for nonlinear time-dependent coupled PDE systems.
//  Solved by NOX solver via Newton, or JFNK with or without preconditioning.
//
//  PDE: Flame propagation (same as tutorial example 17-newton-timedep-flame).
//
//  Domain: Same as in tutorial example 17-newton-timedep-flame.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 2;         // Number of initial uniform mesh refinements.
const int P_INIT = 2;               // Initial polynomial degree of all mesh elements.
const bool JFNK = false;            // true = jacobian-free method,
                                    // false = Newton
const int PRECOND = 1;              // Preconditioning by jacobian (1) (less GMRES iterations, more time to create precond)
                                    // or by approximation of jacobian (2) (less time for precond creation, more GMRES iters).
                                    // in case of jfnk,
                                    // default Ifpack proconditioner in case of Newton.
const double TAU = 0.05;            // Time step.
const bool TRILINOS_OUTPUT = true;  // Display more details about nonlinear and linear solvers.

// Problem parameters.
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;

// Boundary condition types.
BCType bc_types(int marker)
  { return (marker == 1) ? BC_ESSENTIAL : BC_NATURAL; }

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
  { return (ess_bdy_marker == 1) ? 1.0 : 0; }

// Initial conditions.
scalar temp_ic(double x, double y, scalar& dx, scalar& dy)
  { return (x <= x1) ? 1.0 : exp(x1 - x); }

scalar conc_ic(double x, double y, scalar& dx, scalar& dy)
  { return (x <= x1) ? 0.0 : 1.0 - exp(Le*(x1 - x)); }

// Weak forms. 
# include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Create H1 spaces with default shapesets.
  H1Space tspace(&mesh, bc_types, essential_bc_values, P_INIT);
  H1Space cspace(&mesh, bc_types, NULL, P_INIT);
  info("Number of DOF: %d", tspace.get_num_dofs() + cspace.get_num_dofs());

  // Define initial conditions.
  Solution tprev1, cprev1, tprev2, cprev2, titer, citer, tsln, csln;
  tprev1.set_exact(&mesh, temp_ic);  
  cprev1.set_exact(&mesh, conc_ic);
  tprev2.set_exact(&mesh, temp_ic);  
  cprev2.set_exact(&mesh, conc_ic);
  titer.set_exact(&mesh, temp_ic);   
  citer.set_exact(&mesh, conc_ic);
  DXDYFilter omega_dt(omega_dt_fn, &tprev1, &cprev1);
  DXDYFilter omega_dc(omega_dc_fn, &tprev1, &cprev1);

  // Initialize visualization.
  ScalarView rview("Reaction rate", 0, 0, 1600, 460);
  rview.set_min_max_range(0.0,2.0);

  // Initialize weak formulation.
  WeakForm wf(2, JFNK ? true : false);
  if (!JFNK || (JFNK && PRECOND == 1))
  {
    wf.add_jacform(0, 0, callback(jacobian_0_0));
    wf.add_jacform_surf(0, 0, callback(jacobian_0_0_surf));
    wf.add_jacform(1, 1, callback(jacobian_1_1));
    wf.add_jacform(0, 1, callback(jacobian_0_1));
    wf.add_jacform(1, 0, callback(jacobian_1_0));
  }
  else if (PRECOND == 2)
  {
    wf.add_jacform(0, 0, callback(precond_0_0));
    wf.add_jacform(1, 1, callback(precond_1_1));
  }
  wf.add_resform(0, callback(residual_0), H2D_ANY, Tuple<MeshFunction*>(&tprev1, &tprev2));
  wf.add_resform_surf(0, callback(residual_0_surf), 3);
  wf.add_resform(1, callback(residual_1), H2D_ANY, Tuple<MeshFunction*>(&cprev1, &cprev2));

  // Project the functions "titer" and "citer" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solutions on the FE meshes.");

  LinSystem ls(&wf, Tuple<Space*>(&tspace, &cspace));
  ls.project_global(Tuple<MeshFunction*>(&tprev1, &cprev1), 
                    Tuple<Solution*>(&tprev1, &cprev1));

  // Get the coefficient vector.
  scalar *vec = ls.get_solution_vector();

  // Measure the projection time.
  double proj_time = cpu_time.tick().last();

  // Initialize finite element problem.
  H1Shapeset shapeset;
  FeProblem fep(&wf);
  fep.set_spaces(Tuple<Space*>(&tspace, &cspace));
  PrecalcShapeset pss(&shapeset);
  //fep.set_pss(1, &pss);

  // Initialize NOX solver and preconditioner.
  NoxSolver solver(&fep);
  MlPrecond pc("sa");
  if (PRECOND)
  {
    if (JFNK) solver.set_precond(&pc);
    else solver.set_precond("Ifpack");
  }
  if (TRILINOS_OUTPUT)
    solver.set_output_flags(NOX::Utils::Error | NOX::Utils::OuterIteration |
                            NOX::Utils::OuterIterationStatusTest |
                            NOX::Utils::LinearSolverDetails);

  // Time stepping loop:
  double total_time = 0.0;
  cpu_time.tick_reset();
  for (int ts = 1; total_time <= 60.0; ts++)
  {
    info("---- Time step %d, t = %g s", ts, total_time + TAU);

    cpu_time.tick(H2D_SKIP);
    solver.set_init_sln(vec);
    bool solved = solver.solve();
    if (solved)
    {
      vec = solver.get_solution_vector();
      tsln.set_fe_solution(&tspace, &pss, vec);
      csln.set_fe_solution(&cspace, &pss, vec);

      cpu_time.tick();
      info("Number of nonlin iterations: %d (norm of residual: %g)",
          solver.get_num_iters(), solver.get_residual());
      info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)",
          solver.get_num_lin_iters(), solver.get_achieved_tol());

      // Time measurement.
      cpu_time.tick(H2D_SKIP);

      // Visualization.
      DXDYFilter omega_view(omega_fn, &tsln, &csln);
      rview.show(&omega_view);
      cpu_time.tick(H2D_SKIP);
			
      // Skip visualization time.
      cpu_time.tick(H2D_SKIP);

      // Update global time.
      total_time += TAU;

      // Saving solutions for the next time step.
      tprev2.copy(&tprev1);
      cprev2.copy(&cprev1);
      tprev1 = tsln;
      cprev1 = csln;
    }
    else
      error("NOX failed.");

    info("Total running time for time level %d: %g s.", ts, cpu_time.tick().last());
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

