#include "hermes2d.h"
#include "solver_umfpack.h"

using namespace RefinementSelectors;

// This test makes sure that example 11-adapt-system works correctly with MULTI = true.

const int P_INIT = 1;            // Initial polynomial degree of all mesh elements.
const bool MULTI = true;         // MULTI = true  ... use multi-mesh,
                                 // MULTI = false ... use single-mesh.
                                 // Note: In the single mesh option, the meshes are
                                 // forced to be geometrically the same but the
                                 // polynomial degrees can still vary.
const bool SAME_ORDERS = true;   // SAME_ORDERS = true ... when single-mesh is used,
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
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See the Sphinx tutorial (http://hpfem.org/hermes2d/doc/src/tutorial-2.html#adaptive-h-fem-and-hp-fem) for details.
const int MESH_REGULARITY = -1;  // Maximum allowed level of hanging nodes:
                                 // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                 // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                 // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                 // Note that regular meshes are not supported, this is due to
                                 // their notoriously bad performance.
const double CONV_EXP = 1.0;     // Default value is 1.0. This parameter influences the selection of
                                 // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAX_ORDER = 10;        // Maximum allowed element degree
const double ERR_STOP = 5.0;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                 // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;     // Adaptivity process stops when the number of degrees of freedom grows over
                                 // this limit. This is mainly to prevent h-adaptivity to go on forever.

// problem constants
const double E  = 200e9;  // Young modulus for steel: 200 GPa
const double nu = 0.3;    // Poisson ratio
const double f  = 1e3;    // load force: 10^3 N
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));

// boundary markers
const int marker_left = 1;
const int marker_top = 2;

// boundary condition types
int bc_types(int marker)
  { return (marker == marker_left) ? BC_ESSENTIAL : BC_NATURAL; }

// function values for Dirichlet boundary markers
// (if the return value is zero, this can be omitted)
scalar bc_values(int marker, double x, double y)
{
  return 0;
}

// bilinear forms
template<typename Real, typename Scalar>
Scalar bilinear_form_0_0(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (lambda + 2*mu) * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
                      mu * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_0_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return lambda * int_dudy_dvdx<Real, Scalar>(n, wt, u, v) +
             mu * int_dudx_dvdy<Real, Scalar>(n, wt, u, v);
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
  return              mu * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
         (lambda + 2*mu) * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_1(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return -f * int_v<Real, Scalar>(n, wt, v);
}
////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh xmesh, ymesh;
  H2DReader mloader;
  mloader.load("bracket.mesh", &xmesh);

  // initial mesh refinements
  xmesh.refine_element(1);
  xmesh.refine_element(4);

  // create initial mesh for the vertical displacement component,
  // identical to the mesh for the horizontal displacement
  // (bracket.mesh becomes a master mesh)
  ymesh.copy(&xmesh);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset xpss(&shapeset);
  PrecalcShapeset ypss(&shapeset);

  // create the x displacement space
  H1Space xdisp(&xmesh, &shapeset);
  xdisp.set_bc_types(bc_types);
  xdisp.set_bc_values(bc_values);
  xdisp.set_uniform_order(P_INIT);

  // create the y displacement space
  H1Space ydisp(MULTI ? &ymesh : &xmesh, &shapeset);
  ydisp.set_bc_types(bc_types);
  ydisp.set_bc_values(bc_values);
  ydisp.set_uniform_order(P_INIT);

  // enumerate basis functions
  int ndof = assign_dofs(2, &xdisp, &ydisp);

  // initialize the weak formulation
  WeakForm wf(2);
  wf.add_biform(0, 0, callback(bilinear_form_0_0), H2D_SYM);  // note that only one symmetric part is
  wf.add_biform(0, 1, callback(bilinear_form_0_1), H2D_SYM);  // added in the case of symmetric bilinear
  wf.add_biform(1, 1, callback(bilinear_form_1_1), H2D_SYM);  // forms
  wf.add_liform_surf(1, callback(linear_form_surf_1), marker_top);

  // matrix solver
  UmfpackSolver umfpack;

  // create a selector which will select optimal candidate
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, MAX_ORDER, &shapeset);

  // adaptivity loop
  int it = 1;
  bool done = false;
  TimePeriod cpu_time;
  Solution x_sln_coarse, y_sln_coarse;
  Solution x_sln_fine, y_sln_fine;
  do
  {
    info("---- Adaptivity step %d ---------------------------------------------", it); it++;

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // enumerate degrees of freedom
    ndof = assign_dofs(2, &xdisp, &ydisp);

    // solve the coarse mesh problem
    LinSystem ls(&wf, &umfpack);
    ls.set_spaces(2, &xdisp, &ydisp);
    ls.set_pss(2, &xpss, &ypss);
    ls.assemble();
    ls.solve(2, &x_sln_coarse, &y_sln_coarse);

    // time measurement
    cpu_time.tick();

    // report dofs
    info("xdof=%d, ydof=%d\n", xdisp.get_num_dofs(), ydisp.get_num_dofs());

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // solve the fine mesh problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(2, &x_sln_fine, &y_sln_fine);

    // calculate element errors and total error estimate
    H1Adapt hp(Tuple<Space*>(&xdisp, &ydisp));
    hp.set_solutions(Tuple<Solution*>(&x_sln_coarse, &y_sln_coarse), Tuple<Solution*>(&x_sln_fine, &y_sln_fine));
    hp.set_biform(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
    hp.set_biform(0, 1, bilinear_form_0_1<scalar, scalar>, bilinear_form_0_1<Ord, Ord>);
    hp.set_biform(1, 0, bilinear_form_1_0<scalar, scalar>, bilinear_form_1_0<Ord, Ord>);
    hp.set_biform(1, 1, bilinear_form_1_1<scalar, scalar>, bilinear_form_1_1<Ord, Ord>);
    double err_est = hp.calc_error() * 100;

    // time measurement
    cpu_time.tick();

    // report results
    info("Estimate of error: %g%%", err_est);

    // time measurement
    cpu_time.tick(H2D_SKIP);

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY, SAME_ORDERS);
      ndof = assign_dofs(2, &xdisp, &ydisp);
      if (ndof >= NDOF_STOP) done = true;
    }

    // time measurement
    cpu_time.tick();
  }
  while (!done);
  verbose("Total running time: %g s", cpu_time.accumulated());

#define ERROR_SUCCESS       0
#define ERROR_FAILURE      -1
  int ndof_allowed = 820;
  printf("ndof actual = %d\n", ndof);
  printf("ndof allowed = %d\n", ndof_allowed);
  if (ndof <= ndof_allowed) {      // ndofs was 816 at the time this test was created
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
