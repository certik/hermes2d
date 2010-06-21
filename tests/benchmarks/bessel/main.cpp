#include "hermes2d.h"

using namespace RefinementSelectors;

/// This test makes sure that the benchmark "bessel" works correctly.
///
///  Parameters
///  - P_INIT=1
///  - THRESHOLD=0.3
///  - STRATEGY=1
///  - CAND_LIST=HP_ANISO
///  - MESH_REGULARITY=-1
///  - CONV_EXP=1.0
///  - ERR_STOP=0.1
///  - NDOF_STOP=40000
///
///  Results for given parameters
///  - DOFs: 3150

const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.3;     // This is a quantitative parameter of the adapt(...) function and
                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;           // Adaptive strategy:
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
const int MESH_REGULARITY = -1;   // Maximum allowed level of hanging nodes:
                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                  // Note that regular meshes are not supported, this is due to
                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;     // Default value is 1.0. This parameter influences the selection of
                                 // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.1;      // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 40000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// problem constants
const double mu_r   = 1.0;
const double kappa  = 1.0;
const double lambda = 1.0;

// Bessel function of the first kind, order n, defined in bessel.cpp
double jv(double n, double x);

static void exact_sol_val(double x, double y, scalar& e0, scalar& e1)
{
  double t1 = x*x;
  double t2 = y*y;
  double t4 = sqrt(t1+t2);
  double t5 = jv(-1.0/3.0,t4);
  double t6 = 1/t4;
  double t7 = jv(2.0/3.0,t4);
  double t11 = (t5-2.0/3.0*t6*t7)*t6;
  double t12 = atan2(y,x);
  if (t12 < 0) t12 += 2.0*M_PI;
  double t13 = 2.0/3.0*t12;
  double t14 = cos(t13);
  double t17 = sin(t13);
  double t18 = t7*t17;
  double t20 = 1/t1;
  double t23 = 1/(1.0+t2*t20);
  e0 = t11*y*t14-2.0/3.0*t18/x*t23;
  e1 = -t11*x*t14-2.0/3.0*t18*y*t20*t23;
}

static void exact_sol(double x, double y, scalar& e0, scalar& e1, scalar& e1dx, scalar& e0dy)
{
  exact_sol_val(x,y,e0,e1);

  double t1 = x*x;
  double t2 = y*y;
  double t3 = t1+t2;
  double t4 = sqrt(t3);
  double t5 = jv(2.0/3.0,t4);
  double t6 = 1/t4;
  double t7 = jv(-1.0/3.0,t4);
  double t11 = (-t5-t6*t7/3.0)*t6;
  double t14 = 1/t4/t3;
  double t15 = t14*t5;
  double t21 = t7-2.0/3.0*t6*t5;
  double t22 = 1/t3*t21;
  double t27 = atan2(y,x);
  if (t27 < 0) t27 += 2.0*M_PI;
  double t28 = 2.0/3.0*t27;
  double t29 = cos(t28);
  double t32 = t21*t14;
  double t35 = t21*t6;
  double t36 = t35*t29;
  double t39 = sin(t28);
  double t41 = 1/t1;
  double t43 = 1.0+t2*t41;
  double t44 = 1/t43;
  double t47 = 4.0/3.0*t35/x*t39*y*t44;
  double t48 = t5*t29;
  double t49 = t1*t1;
  double t52 = t43*t43;
  double t53 = 1/t52;
  double t57 = t5*t39;
  double t59 = 1/t1/x;
  e1dx =-(t11*x+2.0/3.0*t15*x-2.0/3.0*t22*x)
              *t6*x*t29+t32*t1*t29-t36-t47+4.0/9.0*t48*t2/t49*t53+4.0/3.0*t57*y*t59*t44-4.0/3.0*t57*t2*y/t49/x*t53;
  e0dy = (t11*y+2.0/3.0*t15*y-2.0/3.0*t22*y)*t6*y*t29-t32*t2*t29+t36-t47-4.0/9.0*t48*t41*t53+4.0/3.0*t57*t59*t53*y;
}

// exact solution
scalar2& exact(double x, double y, scalar2& dx, scalar2& dy)
{
  static scalar2 ex;
  exact_sol(x,y, ex[0], ex[1], dx[1], dy[0]);
  return ex;
}

// boundary conditions
BCType bc_types(int marker)
{
  if (marker == 1 || marker == 6)
    return BC_ESSENTIAL; // perfect conductor
  else
    return BC_NATURAL; // impedance
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return 1.0/mu_r * int_curl_e_curl_f<Real, Scalar>(n, wt, u, v) -
         sqr(kappa) * int_e_f<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  cplx ii = cplx(0.0, 1.0);
  return ii * (-kappa) * int_e_tau_f_tau<Real, Scalar>(n, wt, u, v, e);
}

scalar linear_form_surf(int n, double *wt, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
{
  scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    double r = sqrt(e->x[i] * e->x[i] + e->y[i] * e->y[i]);
    double theta = atan2(e->y[i], e->x[i]);
    if (theta < 0) theta += 2.0*M_PI;
    double j13    = jv(-1.0/3.0, r),    j23    = jv(+2.0/3.0, r);
    double cost   = cos(theta),         sint   = sin(theta);
    double cos23t = cos(2.0/3.0*theta), sin23t = sin(2.0/3.0*theta);

    double Etau = e->tx[i] * (cos23t*sint*j13 - 2.0/(3.0*r)*j23*(cos23t*sint + sin23t*cost)) +
                  e->ty[i] * (-cos23t*cost*j13 + 2.0/(3.0*r)*j23*(cos23t*cost - sin23t*sint));

    result += wt[i] * cplx(cos23t*j23, -Etau) * ((v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
  }
  return result;
}
// maximal polynomial order to integrate surface linear form
Ord linear_form_surf_ord(int n, double *wt, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
  {  return Ord(v->val[0].get_max_order());  }


int main(int argc, char* argv[])
{
  // Check input parameters.
  // If true, coarse mesh FE problem is solved in every adaptivity step.
  // If false, projection of the fine mesh solution on the coarse mesh is used. 
  bool SOLVE_ON_COARSE_MESH = false;
  if (argc > 1 && strcmp(argv[1], "-coarse_mesh") == 0)
    SOLVE_ON_COARSE_MESH = true;


  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("lshape3q.mesh", &mesh);

  // Create an Hcurl space with default shapeset.
  HcurlSpace space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf(1);
  wf.add_matrix_form(0, 0, callback(bilinear_form), H2D_SYM);
  wf.add_matrix_form_surf(0, 0, callback(bilinear_form_surf));
  wf.add_vector_form_surf(0, linear_form_surf, linear_form_surf_ord);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu;

  // Initialize refinement selector.
  HcurlProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize the coarse mesh problem.
  LinSystem ls(&wf, &space);

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  Solution sln_coarse, sln_fine;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Assemble and solve the fine mesh problem.
    info("Solving on fine meshes.");
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(&sln_fine);

    // Either solve on coarse mesh or project the fine mesh solution 
    // on the coarse mesh.
    if (SOLVE_ON_COARSE_MESH) {
      info("Solving on coarse mesh.");
      ls.assemble();
      ls.solve(&sln_coarse);
    }
    else {
      info("Projecting fine mesh solution on coarse mesh.");
      int proj_type = 2;    // Hcurl projection.
      ls.project_global(&sln_fine, &sln_coarse, proj_type);
    }

    // Time measurement.
    cpu_time.tick();

    // Calculate element errors and total error estimate.
    info("Calculating error.");
    HcurlAdapt hp(&ls);
    hp.set_solutions(&sln_coarse, &sln_fine); 
    double err_est = hp.calc_error() * 100;
    // double err_est_hcurl = hcurl_error(&sln_coarse, &sln_fine) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est: %g%%", 
      ls.get_num_dofs(), rs.get_num_dofs(), err_est);

    // Add entry to DOF convergence graph.
    graph_dof.add_values(ls.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");

    // Add entry to CPU convergence graph.
    graph_cpu.add_values(cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu.dat");

    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP) done = true;
    else {
      info("Adapting coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      if (ls.get_num_dofs() >= NDOF_STOP) done = true;
    }

    as++;
 }
  while (!done);
  verbose("Total running time: %g s", cpu_time.accumulated());

  int ndof = ls.get_num_dofs();

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  int n_dof_allowed = 3150;
  printf("n_dof_actual = %d\n", ndof);
  printf("n_dof_allowed = %d\n", n_dof_allowed);// ndofs was 2680 at the time this test was created. 
  if (ndof <= n_dof_allowed) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}

