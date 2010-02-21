#include "hermes2d.h"
#include "solver_umfpack.h"

// This example demostrates the adaptive multimesh hp-FEM. The problem
// considered is rooted in linear thermoelasticity: A massive hollow conductor
// is heated by induction and cooled by water running inside. We will approximate the
// x-displacement, y-displacement, and the temperature on individual meshes
// equipped with mutually independent adaptivity mechanisms.
//
// PDE: Linear thermoelasticity
//
// BC: u_1 = u_2 = 0 on Gamma_1
//     du_1/dn = du_2/dn = 0 elsewhere
//     temp = TEMP_INNER on Gamma_4
//     negative heat flux with HEAT_FLUX_OUTER elsewhere
//
// The parameters below can be played with. In particular, compare hp- and
// h-adaptivity via the ADAPT_TYPE option, and compare the multi-mesh vs. single-mesh
// method using the MULTI parameter.

const int P_INIT_TEMP = 2;       // Initial polynomial degrees in temperature mesh.
const int P_INIT_DISP = 2;       // Initial polynomial degrees for displacement meshes.
const int MAX_ORDER = 6;         // Maximum order used during adaptivity.
const bool MULTI = false;         // MULTI = true  ... use multi-mesh,
                                 // MULTI = false ... use single-mesh.
                                 // Note: In the single mesh option, the meshes are
                                 // forced to be geometrically the same but the
                                 // polynomial degrees can still vary.
const bool SAME_ORDERS = false;  // SAME_ORDERS = true ... when single-mesh is used,
                                 // this forces the meshes for all components to be
                                 // identical, including the polynomial degrees of
                                 // corresponding elements. When multi-mesh is used,
                                 // this parameter is ignored.
const double THRESHOLD = 0.3;    // This is a quantitative parameter of the adapt(...) function and
                                 // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;          // Adaptive strategy:
                                 // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                 //   error is processed. If more elements have similar errors, refine
                                 //   all to keep the mesh symmetric.
                                 // STRATEGY = 1 ... refine all elements whose error is larger
                                 //   than THRESHOLD times maximum element error.
                                 // STRATEGY = 2 ... refine all elements whose error is larger
                                 //   than THRESHOLD.
                                 // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int ADAPT_TYPE = 0;        // Type of automatic adaptivity:
                                 // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                 // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                 // ADAPT_TYPE = 2 ... adaptive p-FEM.
const bool ISO_ONLY = false;     // Isotropic refinement flag (concerns quadrilateral elements only).
                                 // ISO_ONLY = false ... anisotropic refinement of quad elements
                                 // is allowed (default),
                                 // ISO_ONLY = true ... only isotropic refinements of quad elements
                                 // are allowed.
const int MESH_REGULARITY = -1;  // Maximum allowed level of hanging nodes:
                                 // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                 // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                 // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                 // Note that regular meshes are not supported, this is due to
                                 // their notoriously bad performance.
const double CONV_EXP = 1.0;     // Default value is 1.0. This parameter influences the selection of 
                                 // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAXIMUM_ORDER = 10;    // Maximum allowed element degree
const double ERR_STOP = 1.0;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                 // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;     // Adaptivity process stops when the number of degrees of freedom grows over
                                 // this limit. This is mainly to prevent h-adaptivity to go on forever.

// problem constants
double HEAT_SRC = 10000.0;       // heat source in the material (caused by induction heating)
double TEMP_INNER = 50;
double HEAT_FLUX_OUTER = -50;
const double E = 2e11;           // steel: E=200 GPa
const double nu = 0.3;
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));
const double l2m = lambda + 2*mu;
const double rho = 8000;
const double g = 9.81;
const double alpha = 13e-6;      // see http://hyperphysics.phy-astr.gsu.edu/hbase/tables/thexp.html

//  Boundary markers:
const int marker_bottom = 1;
const int marker_sides = 2;
const int marker_top = 3;
const int marker_holes = 4;

// boundary markers
int bc_types_x(int marker)
  { return (marker == marker_bottom) ? BC_ESSENTIAL : BC_NATURAL; }

int bc_types_y(int marker)
  { return (marker == marker_bottom) ? BC_ESSENTIAL : BC_NATURAL; }

int bc_types_t(int marker)
  { return (marker == marker_holes) ? BC_ESSENTIAL : BC_NATURAL; }

double bc_values_t(int marker, double x, double y)
  { return TEMP_INNER; }


template<typename Real, typename Scalar>
Scalar bilinear_form_0_0(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return l2m * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
          mu * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_0_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return lambda * int_dudy_dvdx<Real, Scalar>(n, wt, u, v) +
             mu * int_dudx_dvdy<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_0_2(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return - (3*lambda + 2*mu) * alpha * int_dudx_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_0(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return     mu * int_dudy_dvdx<Real, Scalar>(n, wt, u, v) +
         lambda * int_dudx_dvdy<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return  mu * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
         l2m * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_2(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return - (3*lambda + 2*mu) * alpha * int_dudy_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_2_2(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_1(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return -g * rho * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_2(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return HEAT_SRC * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_2(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return HEAT_FLUX_OUTER * int_v<Real, Scalar>(n, wt, v);
}


int main(int argc, char* argv[])
{
  // load the mesh
  Mesh xmesh, ymesh, tmesh;
  H2DReader mloader;
  mloader.load("domain_round_3.mesh", &xmesh); // master mesh
  ymesh.copy(&xmesh);                // ydisp will share master mesh with xdisp
  tmesh.copy(&xmesh);                // temp will share master mesh with xdisp

  // initialize the shapeset and the cache
  H1ShapesetOrtho shapeset;
  PrecalcShapeset xpss(&shapeset);
  PrecalcShapeset ypss(&shapeset);
  PrecalcShapeset tpss(&shapeset);

  // create the x displacement space
  H1Space xdisp(&xmesh, &shapeset);
  xdisp.set_bc_types(bc_types_x);
  xdisp.set_uniform_order(P_INIT_DISP);

  // create the y displacement space
  H1Space ydisp(MULTI ? &ymesh : &xmesh, &shapeset);
  ydisp.set_bc_types(bc_types_y);
  ydisp.set_uniform_order(P_INIT_DISP);

  // create the temperature space
  H1Space temp(MULTI ? &tmesh : &xmesh, &shapeset);
  temp.set_bc_types(bc_types_t);
  temp.set_bc_values(bc_values_t);
  temp.set_uniform_order(P_INIT_TEMP);

  // initialize the weak formulation
  WeakForm wf(3);
  wf.add_biform(0, 0, callback(bilinear_form_0_0));
  wf.add_biform(0, 1, callback(bilinear_form_0_1), SYM);
  wf.add_biform(0, 2, callback(bilinear_form_0_2));
  wf.add_biform(1, 1, callback(bilinear_form_1_1));
  wf.add_biform(1, 2, callback(bilinear_form_1_2));
  wf.add_biform(2, 2, callback(bilinear_form_2_2));
  wf.add_liform(1, callback(linear_form_1));
  wf.add_liform(2, callback(linear_form_2));
  wf.add_liform_surf(2, callback(linear_form_surf_2));

  // visualization
  OrderView xord("X displacement poly degrees", 0, 0, 850, 400);
  OrderView yord("Y displacement poly degrees", 0, 455, 850, 400);
  OrderView tord("Temperature poly degrees", 0, 885, 850, 400);
  ScalarView sview("Von Mises stress [Pa]", 860, 0, 850, 400);
  ScalarView tview("Temperature [deg C]", 860, 455, 850, 400);
  // disable scales
  //xord.show_scale(false); yord.show_scale(false); tord.show_scale(false);

  // matrix solver
  UmfpackSolver solver;

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof, graph_cpu;

  // adaptivity loop
  int it = 1, ndofs;
  bool done = false;
  double cpu = 0.0;
  Solution x_sln_coarse, y_sln_coarse, t_sln_coarse;
  Solution x_sln_fine, y_sln_fine, t_sln_fine;
  do
  {
    info("\n---- Adaptivity step %d ---------------------------------------------\n", it++);

    // time measurement
    begin_time();

    //calculating the number of degrees of freedom
    ndofs = xdisp.assign_dofs(0);
    ndofs += ydisp.assign_dofs(ndofs);
    ndofs += temp.assign_dofs(ndofs);
    printf("xdof=%d, ydof=%d, tdof=%d\n", xdisp.get_num_dofs(), ydisp.get_num_dofs(), temp.get_num_dofs());

    // solve the coarse mesh problem
    LinSystem ls(&wf, &solver);
    ls.set_spaces(3, &xdisp, &ydisp, &temp);
    ls.set_pss(3, &xpss, &ypss, &tpss);
    ls.assemble();
    ls.solve(3, &x_sln_coarse, &y_sln_coarse, &t_sln_coarse);

    // time measurement
    cpu += end_time();

    // view the solution -- this can be slow; for illustration only
    xord.show(&xdisp);
    yord.show(&ydisp);
    tord.show(&temp);
    VonMisesFilter mises(&x_sln_coarse, &y_sln_coarse, mu, lambda);
    sview.set_min_max_range(0, 4e9);
    sview.show(&mises, EPS_HIGH);
    tview.show(&t_sln_coarse, EPS_HIGH);

    // time measurement
    begin_time();

    // solve the fine mesh problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(3, &x_sln_fine, &y_sln_fine, &t_sln_fine);

    // calculate element errors and total error estimate
    H1OrthoHP hp(3, &xdisp, &ydisp, &temp);
    hp.set_biform(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
    hp.set_biform(0, 1, bilinear_form_0_1<scalar, scalar>, bilinear_form_0_1<Ord, Ord>);
    hp.set_biform(0, 2, bilinear_form_0_2<scalar, scalar>, bilinear_form_0_2<Ord, Ord>);
    hp.set_biform(1, 0, bilinear_form_1_0<scalar, scalar>, bilinear_form_1_0<Ord, Ord>);
    hp.set_biform(1, 1, bilinear_form_1_1<scalar, scalar>, bilinear_form_1_1<Ord, Ord>);
    hp.set_biform(1, 2, bilinear_form_1_2<scalar, scalar>, bilinear_form_1_2<Ord, Ord>);
    hp.set_biform(2, 2, bilinear_form_2_2<scalar, scalar>, bilinear_form_2_2<Ord, Ord>);
    double err_est = hp.calc_error_n(3, &x_sln_coarse, &y_sln_coarse, &t_sln_coarse, &x_sln_fine, &y_sln_fine, &t_sln_fine) * 100;

    info("\nEstimate of error: %g%%", err_est);

    // time measurement
    cpu += end_time();

    // add entry to DOF convergence graph
    graph_dof.add_values(x_sln_fine.get_num_dofs() + y_sln_fine.get_num_dofs() + t_sln_fine.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");

    // add entry to CPU convergence graph
    graph_cpu.add_values(cpu, err_est);
    graph_cpu.save("conv_cpu.dat");

    // time measurement
    begin_time();

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY, CONV_EXP, MAXIMUM_ORDER, SAME_ORDERS);
      ndofs = xdisp.assign_dofs();
      ndofs += ydisp.assign_dofs(ndofs);
      if (ndofs >= NDOF_STOP) done = true;
    }

    // time measurement
    cpu += end_time();
  }
  while (!done);
  verbose("Total running time: %g sec", cpu);

  // show the fine solution - this is the final result
  VonMisesFilter stress_fine(&x_sln_fine, &y_sln_fine, mu, lambda);
  sview.set_title("Final solution");
  sview.set_min_max_range(0, 3e4);
  sview.show(&stress_fine);

  // wait for keypress or mouse input
  View::wait("Waiting for all views to be closed.");
  return 0;
};

