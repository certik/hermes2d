#include "hermes2d.h"
#include "solver_umfpack.h"

//  This is another example that allows you to compare h- and hp-adaptivity from the point of view
//  of both CPU time requirements and discrete problem size, look at the quality of the a-posteriori
//  error estimator used by Hermes (exact error is provided), etc. You can also change
//  the parameter MESH_REGULARITY to see the influence of hanging nodes on the adaptive process.
//  The problem is made harder for adaptive algorithms by increasing the parameter SLOPE.
//
//  PDE: -Laplace u = f
//
//  Known exact solution, see functions fn() and fndd()
//
//  Domain: unit square (0, 1)x(0, 1), see the file square.mesh
//
//  BC:  Dirichlet, given by exact solution
//
//  The following parameters can be changed:

const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 4;       // Number of initial uniform mesh refinements.
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
RefinementSelectors::H2DRS_CAND_H_ONLY;         // Type of automatic adaptivity.
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
double EPSILON = 0.01;                // epsilon

// boundary condition types
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Dirichlet boundary condition values
scalar bc_values(int marker, double x, double y)
{
    if (marker == 1)
        return 1;
    else
        return 2-pow(x, 0.1)-pow(y, 0.1);
}

// bilinear form for the Poisson equation
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i=0; i < n; i++) {
        double b_x = 1;
        double b_y = 1;
        result +=
            (b_x * u->val[i] - EPSILON*u->dx[i]) * v->dx[i] +
            (b_y * u->val[i] - EPSILON*u->dy[i]) * v->dy[i];
    }
    return result;
}

// bilinear form for the Poisson equation, stabilization
template<typename Real, typename Scalar>
Scalar bilinear_form_stabilization(int n, double *wt, Func<Real> *u, Func<Real> *v,
        Geom<Real> *e, ExtData<Scalar> *ext)
{
    double h_e = e->element->get_diameter();
    Scalar result = 0;
    for (int i=0; i < n; i++) {
        double b_x = 1;
        double b_y = 1;
        double tau = 1/sqrt(9*pow(4*EPSILON/pow(h_e, 2), 2)+
                pow(2*sqrt(b_x*b_x+b_y*b_y)/h_e, 2));
        //printf("h_e = %f; tau = %f;\n", h_e, tau);
        result +=
            (-b_x * v->dx[i] - b_y * v->dy[i]) * tau *
            (+b_x * u->dx[i] + b_y * u->dy[i]);
    }
    return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_stabilization_order(int n, double *wt, Func<Real> *u, Func<Real> *v,
        Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i=0; i < n; i++) {
        double b_x = 1;
        double b_y = 1;
        double tau = 1;
        result +=
            (-b_x * v->dx[i] - b_y * v->dy[i]) * tau *
            (+b_x * u->dx[i] + b_y * u->dy[i]);
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
        double b_x = 1;
        double b_y = 1;
        // This R makes it nonlinear! So we need to use the Newton method:
        double R = fabs(b_x * u->dx[i] + b_y * u->dy[i]);
        result += s_c * 0.5 * h_e * R *
            (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]) /
                (sqrt(pow(u->dx[i], 2) + pow(u->dy[i], 2)) + 1.e-8);
    }
    return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_shock_capturing_order(int n, double *wt, Func<Real> *u, Func<Real> *v,
        Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i=0; i < n; i++) {
        double b_x = 1;
        double b_y = 1;
        Scalar R = (b_x * u->dx[i] + b_y * u->dy[i]);
        result += R *
            (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]) /
                (sqrt(pow(u->dx[i], 2) + pow(u->dy[i], 2)) + 1.e-8);
    }
    return result;
}

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square_quad.mesh", &mesh);
  // mloader.load("square_tri.mesh", &mesh);
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();
  //mesh.refine_towards_boundary(2, 3);

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
  wf.add_biform(0, 0, bilinear_form_stabilization,
          bilinear_form_stabilization_order);
  //wf.add_biform(0, 0, bilinear_form_shock_capturing,
  //        bilinear_form_shock_capturing_order);

  // visualize solution and mesh
  ScalarView sview("Coarse solution", 0, 0, 500, 400);
  OrderView  oview("Polynomial orders", 505, 0, 500, 400);

  // matrix solver
  UmfpackSolver solver;

  // DOF and CPU convergence graphs
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

    // time measurement
    cpu += end_time();

    // view the solution and mesh
    sview.show(&sln_coarse);
    oview.show(&space);
    sview.wait_for_keypress("wait");
    //break;

    // time measurement
    begin_time();

    // solve the fine mesh problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &sln_fine);

    // calculate error estimate wrt. fine mesh solution
    H1AdaptHP hp(1, &space);
    double err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;
    info("Estimate of error: %g%%", err_est);

    // add entries to DOF convergence graphs
    graph_dof_est.add_values(space.get_num_dofs(), err_est);
    graph_dof_est.save("conv_dof_est.dat");

    // add entries to CPU convergence graphs
    graph_cpu_est.add_values(cpu, err_est);
    graph_cpu_est.save("conv_cpu_est.dat");


    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(THRESHOLD, STRATEGY, &selector, MESH_REGULARITY);
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
