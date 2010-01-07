#include "hermes2d.h"
#include "solver_umfpack.h"

//  This example is a standard nuclear engineering benchmark describing an external-force-driven
//  configuration without fissile materials present, using one-group neutron diffusion approximation.
//
//  Note the way of handling different material parameters. An alternative approach is used in 
//  example "saphir".
//
//  PDE: -div(D(x,y)grad\Phi) + \Sigma_a(x,y)\Phi = Q_{ext}(x,y)
//  where D(x, y) is the diffusion coefficient, \Sigma_a(x,y) the absorption cross-section,
//  and Q_{ext}(x,y) external sources
//  
//  Domain: square (0, L)x(0, L) where L = 30c (see mesh file domain.mesh)
//
//  BC:  Zero Dirichlet for the right and top edges ("vacuum boundary")
//       Zero Neumann for the left and bottom edges ("reflection boundary")
//
//  The following parameters can be changed:

const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 3;       // Number of initial uniform mesh refinements
const double THRESHOLD = 0.6;     // This is a quantitative parameter of the adapt(...) function and
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
const int ADAPT_TYPE = 0;         // Type of automatic adaptivity:
                                  // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                  // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                  // ADAPT_TYPE = 2 ... adaptive p-FEM.
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
const double ERR_STOP = 1e-2;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// Problem parameters
double L = 30;                        // edge of square
double L0 = 0.75*0.5*L;               // end of first water layer
double L1 = 0.5*L;                    // end of second water layer
double L2 = 0.75*L;                   // end of iron layer
double Q_EXT = 1.0;                   // neutron source
double SIGMA_T_WATER = 3.33;
double SIGMA_T_IRON = 1.33;
double C_WATER = 0.994;
double C_IRON = 0.831;

// total cross-section
double sigma_t(double x, double y) {
  if(x < L1 && y < L1) return SIGMA_T_WATER;
  if(x > L2 || y > L2) return SIGMA_T_WATER;
  return SIGMA_T_IRON;
}

// diffusion coefficient
double D(double x, double y) {
  return 1./(3.*sigma_t(x, y));
}

// scattering ratio
double c(double x, double y) {
  if(x < L1 && y < L1) return C_WATER;
  if(x > L2 || y > L2) return C_WATER;
  return C_IRON;
}

// absorption cross section
double sigma_a(double x, double y) {
  double sigma_t_ = sigma_t(x, y);
  double c_ = c(x, y);
  double sigma_s = c_*sigma_t_; // scattering cross-section
  return sigma_t_ - sigma_s;
}

// sources
double q_ext(double x, double y) {
  if(x < L0 && y < L0) return Q_EXT;
  else return 0;
}

/********** Boundary conditions ***********/

// Boundary condition types
int bc_types(int marker)
{
  if (marker == 1) return BC_NATURAL;
  else return BC_ESSENTIAL;
}

// Dirichlet boundary condition values
scalar bc_values(int marker, double x, double y)
{
  return 0.0;
}

// (Volumetric) bilinear form
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    double x = e->x[i];
    double y = e->y[i];
    result += (D(x, y) * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i])
	       + sigma_a(x, y) * u->val[i] * v->val[i]) * wt[i];
  }
  return result;
}

// Integration order for the bilinear form
Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return u->val[0] * v->val[0]; // returning the sum of the degrees of the basis 
                                // and test function
}

// (Volumetrix) linear form (right-hand side)
template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  //return int_F_v<Real, Scalar>(n, wt, q_ext, v, e);
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    double x = e->x[i];
    double y = e->y[i];
    result += q_ext(x, y) * v->val[i] * wt[i];
  }
  return result;
}

// Integration order for the linear form
Ord linear_form_ord(int n, double *wt, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return v->val[0];  // q_ext is piecewise constant, thus 
                     // returning the polynomial degree of the test function;
}

int main(int argc, char* argv[])
{
  // Load the mesh
  Mesh mesh;
  mesh.load("domain.mesh");
  // FIXME: this is temporary, the real benchmark has 10x10 elements
  // and also the values L0, L1 and L2 need to be changed.
  if (INIT_REF_NUM < 3) error("INIT_REF_NUM in this example must be at least 3.");
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // Create finite element space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // Enumerate basis functions
  space.assign_dofs();

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form, bilinear_form_ord, SYM);
  wf.add_liform(0, linear_form, linear_form_ord);

  // Visualize solution and mesh
  ScalarView sview("Coarse solution", 0, 100, 798, 700);
  OrderView  oview("Polynomial orders", 800, 100, 798, 700);

  // Matrix solver
  UmfpackSolver solver;

  // Convergence graph wrt. the number of degrees of freedom
  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Error Convergence for the Iron-Water Problem", "Degrees of Freedom", "Error Estimate [%]");
  graph.add_row("error estimate", "k", "-", "o");

  // Convergence graph wrt. CPU time
  GnuplotGraph graph_cpu;
  graph_cpu.set_captions("Error Convergence for the Iron-Water Problem", "CPU Time", "Error Estimate [%]");
  graph_cpu.add_row("error estimate", "k", "-", "o");
  graph_cpu.set_log_y();

  // Adaptivity loop
  int it = 1, ndofs;
  bool done = false;
  double cpu = 0.0;
  Solution sln_coarse, sln_fine;
  do
    {
    info("\n---- Adaptivity step %d ---------------------------------------------\n", it++);

    // Time measurement
    begin_time();

    // Solve the coarse mesh problem
    LinSystem ls(&wf, &solver);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln_coarse);

    // Time measurement
    cpu += end_time();

    // View the solution and mesh
    sview.show(&sln_coarse);
    oview.show(&space);
    //oview.wait_for_keypress();

    // Time measurement
    begin_time();

    //break;

    // Solve the fine mesh problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &sln_fine);

    // Calculate error estimate wrt. fine mesh solution
    H1OrthoHP hp(1, &space);
    double err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;
    info("Estimate of error: %g%%", err_est);

    // Add entry to DOF convergence graph
    graph.add_values(0, space.get_num_dofs(), err_est);
    graph.save("conv_dof.gp");

    // Add entry to CPU convergence graph
    graph_cpu.add_values(0, cpu, err_est);
    graph_cpu.save("conv_cpu.gp");

    // If err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY);
      ndofs = space.assign_dofs();
      if (ndofs >= NDOF_STOP) done = true;
    }

    // Time measurement
    cpu += end_time();
  }
  while (done == false);
  verbose("Total running time: %g sec", cpu);

  // Show the fine solution - this is the final result
  sview.set_title("Final solution");
  sview.show(&sln_fine);

  // Wait for keyboard or mouse input
  View::wait("Waiting for keyboard or mouse input.");
  return 0;
}
