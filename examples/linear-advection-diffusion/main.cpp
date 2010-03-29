#include "hermes2d.h"
#include "solver_umfpack.h"


//  This example solves a linear advection diffusion problem (using optional 
//  variational multiscale stabilization.
//
//  PDE: div(bu - \epsilon \nabla u) = 0 where b = (b1, b2) is a constant vector
//
//  Domain: Square (0, 1)x(0, 1)
//
//  BC:  Dirichlet, see the function scalar bc_values() below.
//
//  The following parameters can be changed:

const int P_INIT = 2;                   // Initial polynomial degree of all mesh elements.
const bool STABILIZATION_ON = false;    // Stabilization on/off.
const bool SHOCK_CAPTURING_ON = false;  // Shock capturing on/off.
const int INIT_REF_NUM = 1;       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 5;   // Number of initial uniform mesh refinements in the boundary layer region.
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
const RefinementSelectors::AllowedCandidates ADAPT_TYPE =
RefinementSelectors::H2DRS_CAND_HP;         // Type of automatic adaptivity.
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
const double CONV_EXP = 1.0;      // Default value is 1.0. This parameter influences the selection of
                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.1;      // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// problem constants
const double EPSILON = 0.01;      // diffusivity           
const double B1 = 1., B2 = 1.;    // advection direction, div(B) = 0

// boundary condition types
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Dirichlet boundary condition values
scalar bc_values(int marker, double x, double y)
{
    if (marker == 1) return 1;
    else return 2 - pow(x, 0.1) - pow(y, 0.1);
}

// bilinear form
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  //result += EPSILON * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
  for (int i=0; i < n; i++) {
    // int_grad_u_grad_v:
    result += wt[i] * EPSILON * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
    // VERSION A:
    //result += wt[i] * (B1 * u->dx[i] * v->val[i] + B2 * u->dy[i] * v->val[i]);
    // VERSION B:
    result += - wt[i] * (B1 * u->val[i] * v->dx[i] + B2 * u->val[i] * v->dy[i]);
  }
  return result;
}

// bilinear form for the stabilization
template<typename Real, typename Scalar>
Scalar bilinear_form_stabilization(int n, double *wt, Func<Real> *u, Func<Real> *v,
        Geom<Real> *e, ExtData<Scalar> *ext)
{
  double h_e = e->element->get_diameter();
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    double b_norm = sqrt(B1*B1 + B2*B2);
    double tau = 1. / sqrt(9*pow(4*EPSILON/pow(h_e, 2), 2) + pow(2*b_norm/h_e, 2));
    // FIXME: we are missing the Laplace operator in both the direct and adjoint operator,
    // waiting for second derivatives
    result += wt[i]*(B1 * v->dx[i] + B2 * v->dy[i]) * tau * (B1 * u->dx[i] + B2 * u->dy[i]);
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_shock_capturing(int n, double *wt, Func<Real> *u, Func<Real> *v,
        Geom<Real> *e, ExtData<Scalar> *ext)
{
  double h_e = e->element->get_diameter();
  double s_c = 0.9;
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    // This R makes it nonlinear! So we need to use the Newton method:
    double R = fabs(B1 * u->dx[i] + B2 * u->dy[i]);
    result += wt[i] * s_c * 0.5 * h_e * R *
              (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]) /
              (sqrt(pow(u->dx[i], 2) + pow(u->dy[i], 2)) + 1.e-8);
  }
  return result;
}

Ord bilinear_form_shock_capturing_order(int n, double *wt, Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext)
{
  return (B1 * u->dx[0] + B2 * u->dy[0]) * (u->dx[0]*v->dx[0] + u->dy[0]*v->dy[0]) /
         (sqrt(pow(u->dx[0], 2) + pow(u->dy[0], 2)) + 1.e-8);
}

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square_quad.mesh", &mesh);
  // mloader.load("square_tri.mesh", &mesh);
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(2, INIT_REF_NUM_BDY);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create finite element space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // enumerate basis functions
  space.assign_dofs();

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form));
  if (STABILIZATION_ON == true) {
    wf.add_biform(0, 0, callback(bilinear_form_stabilization));
  }
  if (SHOCK_CAPTURING_ON == true) {
    wf.add_biform(0, 0, bilinear_form_shock_capturing, bilinear_form_shock_capturing_order);
  }

  // visualize solution and mesh
  OrderView  oview("Coarse mesh", 0, 0, 500, 400);
  ScalarView sview("Coarse mesh solution", 510, 0, 500, 400);
  ScalarView sview2("Fine mesh solution", 1020, 0, 500, 400);

  // matrix solver
  UmfpackSolver solver;

  // DOF convergence graph
  SimpleGraph graph_dof_est, graph_cpu_est;

  // prepare selector
  RefinementSelectors::H1NonUniformHP selector(ISO_ONLY, ADAPT_TYPE, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);

  // adaptivity loop
  int it = 1, ndofs;
  bool done = false;
  double cpu = 0.0;
  Solution sln_coarse, sln_fine;
  do
  {
    info("\n---- Adaptivity step %d ---------------------------------------------\n", it++);

    // time measurement
    begin_time();

    // solve the coarse mesh problem
    LinSystem ls(&wf, &solver);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln_coarse);

    // solve the fine mesh problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &sln_fine);

    // time measurement
    cpu += end_time();

    // show fine mesh solution
    // view the solution and mesh
    oview.show(&space);
    sview.show(&sln_coarse);
    sview2.show(&sln_fine);
    //sview.wait_for_keypress();
    //break;
    
    // time measurement
    begin_time();

    // calculate error estimate wrt. fine mesh solution
    H1AdaptHP hp(1, &space);
    double err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;
    info("Estimate of error: %g%%", err_est);

    // add entries to DOF and CPU convergence graphs
    graph_dof_est.add_values(space.get_num_dofs(), err_est);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu, err_est);
    graph_cpu_est.save("conv_cpu_est.dat");

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(THRESHOLD, STRATEGY, &selector, MESH_REGULARITY);
      //mesh.refine_all_elements();
      //space.set_uniform_order(P_INIT);
      ndofs = space.assign_dofs();
      if (ndofs >= NDOF_STOP) done = true;
    }

    // time measurement
    cpu += end_time();
  }
  while (done == false);
  verbose("Total running time: %g sec", cpu);

  // show the fine solution - this is the final result
  sview.set_title("Final solution");
  sview.show(&sln_fine);

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}
