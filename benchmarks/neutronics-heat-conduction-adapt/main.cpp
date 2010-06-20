#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include <cmath>
#include <iostream>

using namespace RefinementSelectors;

// Neutronics/heat conduction test case (adaptive).
//
// Author: Damien Lebrun-Grandie (Texas A&M University).
//
// PDE:
//
// 1/v d/dt phi = div(D grad phi) + nu Sigma_f phi_1 + q
//
// rho c_p d/dt T = div(k grad T) + kappa Sigma_f phi + qT
//
// Domain: rectangle (Lx, Ly).
//
// BC: homogeneous Dirichlet.
//
// The following parameters can be changed:

const bool SOLVE_ON_COARSE_MESH = true;   // true... Newton is done on coarse mesh in every adaptivity step.
                                           // false...Newton is done on coarse mesh only once, then projection
                                           // of the fine mesh solution to coarse mesh is used.
const int INIT_GLOB_REF_NUM = 2;           // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;            // Number of initial refinements towards boundary.
const int P_INIT = 2;                      // Initial polynomial degree of all mesh elements

// Time-stepping:
const double TAU = 0.1;                    // Time step.
const double T_FINAL = 10.0;               // Time interval length.

// Adaptivity:
const int UNREF_FREQ = 1;                  // Every UNREF_FREQ time step the mesh is unrefined.
const double THRESHOLD = 0.3;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                    // Adaptive strategy:
                                           // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                           //   error is processed. If more elements have similar errors, refine
                                           //   all to keep the mesh symmetric.
                                           // STRATEGY = 1 ... refine all elements whose error is larger
                                           //   than THRESHOLD times maximum element error.
                                           // STRATEGY = 2 ... refine all elements whose error is larger
                                           //   than THRESHOLD.
                                           // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;   // Predefined list of element refinement candidates. Possible values are
                                           // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                           // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                           // See User Documentation for details.
const int MESH_REGULARITY = -1;            // Maximum allowed level of hanging nodes:
                                           // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                           // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                           // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                           // Note that regular meshes are not supported, this is due to
                                           // their notoriously bad performance.
const double CONV_EXP = 1.0;               // Default value is 1.0. This parameter influences the selection of
                                           // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAX_P = 6;                       // Maximum polynomial order allowed in hp-adaptivity
                                           // had to be limited due to complicated integrals.
const double ERR_STOP = 0.01;              // Stopping criterion for hp-adaptivity
                                           // (relative error between reference and coarse solution in percent).
const int NDOF_STOP = 10000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.


// Newton's method:
const double NEWTON_TOL_COARSE = 1.0e-2;   // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 5.0e-2;     // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;            // Maximum allowed number of Newton iterations.

// Problem parameters.
const double CT = 1.0;
const double CF = 1.0;
const double rT = 1.0;
const double rF = 0.25;
const double LX = 100.0;          // Domain sizes in the x and y dimensions.
const double LY = 100.0;
const double invvel = 2.0e-4;     // Inverse of neutron velocity.
const double xsdiff = 1.268;      // Diffusion coefficient.
const double Tref = 0.0;          // Temperature at boundary.

const double nu = 2.41;           // Number of neutrons emitted per fission event.
const double xsfiss = 0.00191244; // Fission cross section.
const double kappa = 1.0e-6;
const double rho = 1.0;           // Density.
const double cp = 1.0;            // Heat capacity.

const double PI = acos(-1.0);
const double normalization_const = 1.0;

const double energy_per_fission = kappa * xsfiss;

// Miscellaneous:
double TIME = 0.0;                // Current time.

// Thermal conductivity depends on temperature
const  double k0 = 3.0e-3;
const  double k1 = 2.0e-4;
template<typename Real>
Real k(Real T) {
  return k0 + k1 * (T - Tref);
}

// Derivative of the thermal conductivity
template<typename Real>
Real dk_dT(Real T) {
  return k1;
}

// Removal cross section depends on temperature
const double xsa_ref = 0.0349778;
const double doppler_coeff = 1.0e-5;
template<typename Real>
Real xsrem(Real T) {
  return xsa_ref + doppler_coeff * (sqrt(T+1.0e-10) - sqrt(Tref));
  //return xsa_ref + doppler_coeff * (sqrt(T) - sqrt(Tref));
}

// Derivative of the removal cross section with respect to temperature
template<typename Real>
Real dxsrem_dT(Real T) {
  return doppler_coeff / (2*sqrt(T+1.0e-10));
  //return doppler_coeff / (2*sqrt(T));
}

// Heat source.
template<typename Real>
Real qT(Real x, Real y) {
  return
rho*cp*CT*(1.0-pow(tanh(rT*TIME),2.0))*rT*sin(x/LX*PI)*sin(y/LY*PI)-k1*CT*CT*pow(1.0+tanh(rT*TIME),2.0)*pow(cos(x/LX*PI),2.0)/(LX*LX)*PI*PI*pow(sin(y/LY*PI),2.0)+(k0+k1*(CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI)-Tref))*CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)/(LX*LX)*PI*PI*sin(y/LY*PI)-k1*CT*CT*pow(1.0+tanh(rT*TIME),2.0)*pow(sin(x/LX*PI),2.0)*pow(cos(y/LY*PI),2.0)/(LY*LY)*PI*PI+(k0+k1*(CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI)-Tref))*CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI)/(LY*LY)*PI*PI-normalization_const*energy_per_fission*xsfiss*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY;
}

// Extraneous neutron source.
template<typename Real>
Real q(Real x, Real y) {
  return 
invvel*CF*rF*exp(rF*TIME)*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY-xsdiff*(-CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)/(LX*LX*LX)*PI*PI*sin(y/LY*PI)*x*y/LY+2.0*CF*(1.0+exp(rF*TIME))*cos(x/LX*PI)/(LX*LX)*PI*sin(y/LY*PI)*y/LY)-xsdiff*(-CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)/(LY*LY*LY)*PI*PI*x/LX*y+2.0*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*cos(y/LY*PI)/(LY*LY)*PI*x/LX)+(xsa_ref+doppler_coeff*(sqrt(CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI))-sqrt(Tref)))*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY-nu*xsfiss*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY;
}

// Boundary condition types.
BCType bc_types_T(int marker)
{
  return BC_ESSENTIAL;
}

BCType bc_types_phi(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values_T(int ess_bdy_marker, double x, double y)
{
  return Tref;
}
 
scalar essential_bc_values_phi(int ess_bdy_marker, double x, double y)
{
  return 0.0;
}

// Weak forms.
# include "forms.cpp"

// Exact solutions.
#include "exact_solution.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh basemesh, mesh_T, mesh_phi;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh);

  // Perform initial mesh refinements.
  for (int i=0; i < INIT_GLOB_REF_NUM; i++)
    basemesh.refine_all_elements();
  basemesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);
  mesh_T.copy(&basemesh);
  mesh_phi.copy(&basemesh);

  // Create H1 spaces with default shapesets.
  H1Space space_T(&mesh_T, bc_types_T, essential_bc_values_T, P_INIT);
  H1Space space_phi(&mesh_phi, bc_types_phi, essential_bc_values_phi, P_INIT);

  // Solutions for the Newton's iteration and time stepping.
  Solution T_prev_newton, T_prev_time, phi_prev_newton, phi_prev_time;

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, jac_TT, jac_TT_ord, H2D_UNSYM, H2D_ANY, &T_prev_newton);
  wf.add_matrix_form(0, 1, jac_Tphi, jac_Tphi_ord, H2D_UNSYM, H2D_ANY, 0);
  wf.add_vector_form(0, res_T, res_T_ord, H2D_ANY, 
                Tuple<MeshFunction*>(&T_prev_newton, &T_prev_time, &phi_prev_newton));
  wf.add_matrix_form(1, 0, jac_phiT, jac_phiT_ord, H2D_UNSYM, H2D_ANY, 
                Tuple<MeshFunction*>(&phi_prev_newton, &T_prev_newton));
  wf.add_matrix_form(1, 1, jac_phiphi, jac_phiphi_ord, H2D_UNSYM, H2D_ANY, &T_prev_newton);
  wf.add_vector_form(1, res_phi, res_phi_ord, H2D_ANY, 
                Tuple<MeshFunction*>(&phi_prev_newton, &phi_prev_time, &T_prev_newton));

  // Initialize solution views.
  ScalarView view_T("", 360, 0, 350, 250);
  view_T.fix_scale_width(80);
  ScalarView view_T_exact("", 0, 0, 350, 250);
  view_T_exact.fix_scale_width(80);
  ScalarView view_phi("", 360, 300, 350, 250);
  view_phi.fix_scale_width(80);
  ScalarView view_phi_exact("", 0, 300, 350, 250);
  view_phi_exact.fix_scale_width(80);

  // Initialize mesh views.
  OrderView ordview_T_coarse("", 720, 0, 350, 250);
  ordview_T_coarse.fix_scale_width(80);
  OrderView ordview_T_fine("", 1080, 0, 350, 250);
  ordview_T_fine.fix_scale_width(80);
  OrderView ordview_phi_coarse("", 720, 300, 350, 250);
  ordview_phi_coarse.fix_scale_width(80);
  OrderView ordview_phi_fine("", 1080, 300, 350, 250);
  ordview_phi_fine.fix_scale_width(80);

  // Initialize the nonlinear system.
  NonlinSystem nls(&wf, Tuple<Space*>(&space_T, &space_phi));

  // Project initial conditions on FE spaces to obtain initial 
  // vector for the Newton's method.
  info("Projecting initial conditions to obtain initial vector for the Newton's method.");
  T_prev_time.set_exact(&mesh_T, T_exact);
  phi_prev_time.set_exact(&mesh_phi, phi_exact);
  nls.project_global(Tuple<MeshFunction*>(&T_prev_time, &phi_prev_time), 
                     Tuple<Solution*>(&T_prev_newton, &phi_prev_newton));

  // Newton's loop on the coarse mesh.
  info("Solving on coarse meshes.");
  bool verbose = true; // Default is false.
  if (!nls.solve_newton(Tuple<Solution*>(&T_prev_newton, &phi_prev_newton), 
                        NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
    error("Newton's method did not converge.");

  // Store the result in T_coarse, phi_coarse.
  Solution T_coarse, phi_coarse;
  T_coarse.copy(&T_prev_newton);
  phi_coarse.copy(&phi_prev_newton);

  // Exact solutions for error evaluation.
  ExactSolution T_solution(&mesh_T, T_exact);
  ExactSolution phi_solution(&mesh_phi, phi_exact);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Time stepping loop:
  Solution T_fine, phi_fine;
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= nstep; ts++)
  {
    // Update global time.
    TIME = ts*TAU;

    // Update time-dependent exact solutions.
    T_solution.update(&mesh_T, T_exact);
    phi_solution.update(&mesh_phi, phi_exact);

    // Show exact solution.
    char title[100];
    sprintf(title, "T (exact), t = %g s", TIME);
    view_T_exact.set_title(title);
    view_T_exact.show_mesh(false);
    view_T_exact.show(&T_solution);
    sprintf(title, "phi (exact), t = %g s", TIME);
    view_phi_exact.set_title(title);
    view_phi_exact.show_mesh(false);
    view_phi_exact.show(&phi_solution);

    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("---- Time step %d - prior to adaptivity:", ts);
      info("Global mesh derefinement.");
      mesh_T.copy(&basemesh);
      mesh_phi.copy(&basemesh);
      space_T.set_uniform_order(P_INIT);
      space_phi.set_uniform_order(P_INIT);

      // Project fine mesh solutions on globally derefined meshes.
      if (SOLVE_ON_COARSE_MESH) 
        info("Projecting fine mesh solution to obtain initial vector on globally derefined mesh.");
      else 
        info("Projecting fine mesh solution on globally derefined mesh for error calculation.");
      nls.project_global(Tuple<MeshFunction*>(&T_fine, &phi_fine), 
                         Tuple<Solution*>(&T_prev_newton, &phi_prev_newton));

      if (SOLVE_ON_COARSE_MESH) {
        // Newton's loop on the globally derefined mesh.
        info("Solving on globally derefined meshes.", ts);
        if (!nls.solve_newton(Tuple<Solution*>(&T_prev_newton, &phi_prev_newton), 
                              NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
          error("Newton's method did not converge.");
      }

      // Store the result in T_coarse, phi_coarse.
      T_coarse.copy(&T_prev_newton);
      phi_coarse.copy(&phi_prev_newton);
    }

    // Adaptivity loop:
    bool done = false;
    double err_est;
    int as = 1;
    do {
      info("---- Time step %d, adaptivity step %d:", ts, as);

      // Initialize reference nonlinear system.
      RefSystem rnls(&nls);

      // Set initial condition for the Newton's method on the fine mesh.
      if (as == 1) {
        info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
        rnls.project_global(Tuple<MeshFunction*>(&T_coarse, &phi_coarse), 
                            Tuple<Solution*>(&T_prev_newton, &phi_prev_newton));
      }
      else {
        info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
        rnls.project_global(Tuple<MeshFunction*>(&T_fine, &phi_fine), 
                            Tuple<Solution*>(&T_prev_newton, &phi_prev_newton));
      }

      // Newton's loop on the fine meshes.
      info("Solving on fine meshes.");
      if (!rnls.solve_newton(Tuple<Solution*>(&T_prev_newton, &phi_prev_newton), 
                             NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose))
        error("Newton's method did not converge.");

      // Store the result in T_fine, phi_fine.
      T_fine.copy(&T_prev_newton);
      phi_fine.copy(&phi_prev_newton);

      // Visualize intermediate solutions and mesh during adaptivity.
      sprintf(title, "T (fine mesh), t = %g s, adapt step %d", TIME, as);
      view_T.set_title(title);
      view_T.show(&T_fine);
      sprintf(title, "phi (fine mesh), t = %g s, adapt step %d", TIME, as);
      view_phi.set_title(title);
      view_phi.show(&phi_fine);
      sprintf(title, "T mesh (coarse), t = %g, adapt step %d", TIME, as);
      ordview_T_coarse.set_title(title);
      ordview_T_coarse.show(&space_T);
      sprintf(title, "T mesh (fine), t = %g, adapt step %d", TIME, as);
      ordview_T_fine.set_title(title);
      ordview_T_fine.show(rnls.get_space(0));
      sprintf(title, "phi mesh (coarse), t = %g, adapt step %d", TIME, as);
      ordview_phi_coarse.set_title(title);
      ordview_phi_coarse.show(&space_phi);
      sprintf(title, "phi mesh (fine), t = %g, adapt step %d", TIME, as);
      ordview_phi_fine.set_title(title);
      ordview_phi_fine.show(rnls.get_space(1));

      // Calculate error estimates and exact errors.
      info("Calculating errors.");
      double T_err_est = h1_error(&T_coarse, &T_fine) * 100;
      double phi_err_est = h1_error(&phi_coarse, &phi_fine) * 100;
      double T_err_exact = h1_error(&T_coarse, &T_solution) * 100;
      double phi_err_exact = h1_error(&phi_coarse, &phi_solution) * 100;
      info("T: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%, err_exact: %g %%", 
	   nls.get_num_dofs(0), rnls.get_num_dofs(0), T_err_est, T_err_exact);
      info("phi: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%, err_exact: %g %%", 
	   nls.get_num_dofs(1), rnls.get_num_dofs(1), phi_err_est, phi_err_exact);
 
      // Calculate element errors and total error estimate for adaptivity.
      H1Adapt hp(&nls);
      hp.set_solutions(Tuple<Solution*>(&T_coarse, &phi_coarse), 
                       Tuple<Solution*>(&T_fine, &phi_fine));
      err_est = hp.calc_error() * 100;

      // If err_est too large, adapt the mesh.
      if (err_est < ERR_STOP) done = true;
      else {
        info("Adapting the coarse meshes.");
        done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (nls.get_num_dofs() >= NDOF_STOP) {
          done = true; 
          break;
        }

        // Project the fine mesh solutions on the new coarse meshes.
        if (SOLVE_ON_COARSE_MESH) 
          info("Projecting fine mesh solution to obtain initial vector on new coarse mesh.");
        else 
          info("Projecting fine mesh solution on coarse mesh for error calculation.");
        nls.project_global(Tuple<MeshFunction*>(&T_fine, &phi_fine), 
                           Tuple<Solution*>(&T_prev_newton, &phi_prev_newton));

        if (SOLVE_ON_COARSE_MESH) {
          // Newton's loop on the coarse meshes.
          info("Solving on coarse meshes.");
          if (!nls.solve_newton(Tuple<Solution*>(&T_prev_newton, &phi_prev_newton), 
                                NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
            error("Newton's method did not converge.");
        }

        // Store the results in T_coarse, phi_coarse.
        T_coarse.copy(&T_prev_newton);
        phi_coarse.copy(&phi_prev_newton);

        as++;
      }
    }
    while (!done);

    // Copy new time level solutions into T_prev_time, phi_prev_time.
    T_prev_time.copy(&T_fine);
    phi_prev_time.copy(&phi_fine);

    // Compute exact error.
    double T_error = h1_error(&T_prev_time, &T_solution) * 100;
    double phi_error = h1_error(&phi_prev_time, &phi_solution) * 100;
    info("Exact solution error for T (H1 norm): %g %%", T_error);
    info("Exact solution error for phi (H1 norm): %g %%", phi_error);
  }


  // Wait for all views to be closed.
  View::wait();
  return 0;
}
