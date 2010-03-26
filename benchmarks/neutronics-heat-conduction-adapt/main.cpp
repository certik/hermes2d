#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"
#include <cmath>
#include <iostream>

using namespace RefinementSelectors;

// Neutronics/Heat Conduction Test Case
//
// -->PDEs:
//
// 1/v d/dt phi = div(D grad phi) + nu Sigma_f phi_1 + q
//
// rho c_p d/dt T = div(k grad T) + kappa Sigma_f phi + qT
//
// -->Domain: rectangle Lx Ly
//
// -->BC: homogeneous Dirichlet
//



const int PROJ_TYPE = 1;                 // For the projection of the initial condition 
                                         // on the initial mesh: 1 = H1 projection, 0 = L2 projection

// Initial refinements
const int INIT_GLOB_REF_NUM = 3;         // Number of initial uniform mesh refinements
const int INIT_BDY_REF_NUM = 0;          // Number of initial refinements towards boundary


// Time-stepping
double TIME = 0.0;
const double TAU = 0.1;                  // Time step
const double T_FINAL = 10.0;             // Time interval length

const int P_INIT = 1;                    // Initial polynomial degree of all mesh elements

// Adaptivity
const int UNREF_FREQ = 1;                // Every UNREF_FREQ time step the mesh is unrefined
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
                                         // See the Sphinx tutorial (http://hpfem.org/hermes2d/doc/src/tutorial-2.html#adaptive-h-fem-and-hp-fem) for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAX_P = 6;                     // Maximum polynomial order allowed in hp-adaptivity
                                         // had to be limited due to complicated integrals
const double ERR_STOP = 0.01;            // Stopping criterion for hp-adaptivity
                                         // (relative error between reference and coarse solution in percent)
const int NDOF_STOP = 10000;             // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.


// Newton's method
const double NEWTON_TOL_COARSE = 1.0e-2; // Stopping criterion for Newton on coarse mesh
const double NEWTON_TOL_FINE = 5.0e-2;   // Stopping criterion for Newton on fine mesh
const int NEWTON_MAX_ITER = 20;          // Maximum allowed number of Newton iterations



/*********************************************************************************************************/
/********** Problem parameters ***********/


const double CT = 1.0;
const double CF = 1.0;
const double rT = 1.0;
const double rF = 0.25;
const double LX = 100.0;          // Domain sizes in the x and y dimensions
const double LY = 100.0;
const double invvel = 2.0e-4;     // Inverse of neutron velocity
const double xsdiff = 1.268;      // Diffusion coefficient
const double Tref = 0.0;         // Temperature at boundary

const double nu = 2.41;           // Number of neutrons emitted per fission event
const double xsfiss = 0.00191244; // Fission cross section
const double kappa = 1.0e-6;
const double rho = 1.0;           // Density
const double cp = 1.0;            // Heat capacity

const double PI = acos(-1.0);
const double normalization_const = 1.0;

const double energy_per_fission = kappa * xsfiss;

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

// Heat source
template<typename Real>
Real qT(Real x, Real y) {
  //  std::cout<<"entering qt"<<std::endl;
  return
rho*cp*CT*(1.0-pow(tanh(rT*TIME),2.0))*rT*sin(x/LX*PI)*sin(y/LY*PI)-k1*CT*CT*pow(1.0+tanh(rT*TIME),2.0)*pow(cos(x/LX*PI),2.0)/(LX*LX)*PI*PI*pow(sin(y/LY*PI),2.0)+(k0+k1*(CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI)-Tref))*CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)/(LX*LX)*PI*PI*sin(y/LY*PI)-k1*CT*CT*pow(1.0+tanh(rT*TIME),2.0)*pow(sin(x/LX*PI),2.0)*pow(cos(y/LY*PI),2.0)/(LY*LY)*PI*PI+(k0+k1*(CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI)-Tref))*CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI)/(LY*LY)*PI*PI-normalization_const*energy_per_fission*xsfiss*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY;
}

// Extraneous neutron source
template<typename Real>
Real q(Real x, Real y) {
  //  std::cout<<"entering q"<<std::endl;
  return 
invvel*CF*rF*exp(rF*TIME)*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY-xsdiff*(-CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)/(LX*LX*LX)*PI*PI*sin(y/LY*PI)*x*y/LY+2.0*CF*(1.0+exp(rF*TIME))*cos(x/LX*PI)/(LX*LX)*PI*sin(y/LY*PI)*y/LY)-xsdiff*(-CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)/(LY*LY*LY)*PI*PI*x/LX*y+2.0*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*cos(y/LY*PI)/(LY*LY)*PI*x/LX)+(xsa_ref+doppler_coeff*(sqrt(CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI))-sqrt(Tref)))*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY-nu*xsfiss*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY;
}

/*********************************************************************************************************/
/********** Boundary conditions ***********/


// Boundary condition types (essential = Dirichlet)
int bc_types_T(int marker)
{
  return BC_ESSENTIAL;
}
 

// Dirichlet boundary condition values
scalar bc_values_T(int marker, double x, double y)
{
  return Tref;
}


int bc_types_phi(int marker)
{
  return BC_ESSENTIAL;
}
 
scalar bc_values_phi(int marker, double x, double y)
{
  return 0.0;
}

/*********************************************************************************************************/
/********** Jacobian and residual  ***********/


# include "forms.cpp"

/*********************************************************************************************************/
/********** Initial conditions ***********/


scalar T_ic(double x, double y, double& dx, double& dy)
{ 
  dx = 1.0;  dy = 0.0;  return x; 
}

scalar phi_ic(double x, double y, double& dx, double& dy)
{
  dx = 0.0;  dy = 1.0;  return y;
}

/*********************************************************************************************************/
/********** Exact solutions ***********/


scalar phi_exact(double x, double y, double& dx, double& dy)
{
  dx = CF*(1+exp(rF*TIME))*sin(PI*x/LX)*sin(PI*y/LY)/LX*y/LY
       + CF*(1+exp(rF*TIME))*PI/LX*cos(PI*x/LX)*sin(PI*y/LY)*x/LX*y/LY;
  dy = CF*(1+exp(rF*TIME))*sin(PI*x/LX)*sin(PI*y/LY)*x/LX/LY
       + CF*(1+exp(rF*TIME))*sin(PI*x/LX)*PI/LY*cos(PI*y/LY)*x/LX*y/LY;
  return CF*(1+exp(rF*TIME))*sin(PI*x/LX)*sin(PI*y/LY)*x/LX*y/LY;
}

scalar T_exact(double x, double y, double& dx, double& dy)
{
  dx = CT*(1+tanh(rT*TIME))*PI/LX*cos(PI*x/LX)*sin(PI*y/LY);
  dy = CT*(1+tanh(rT*TIME))*sin(PI*x/LX)*PI/LY*cos(PI*y/LY);
  return CT*(1+tanh(rT*TIME))*sin(PI*x/LX)*sin(PI*y/LY);
}

/*********************************************************************************************************/
/********** This is the main ***********/


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh basemesh, mesh_T, mesh_phi;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh);

  // perform initial refinements
  for (int i=0; i < INIT_GLOB_REF_NUM; i++)
    basemesh.refine_all_elements();
  basemesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);
  mesh_T.copy(&basemesh);
  mesh_phi.copy(&basemesh);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create H1 spaces
  H1Space space_T(&mesh_T, &shapeset);
  space_T.set_bc_types(bc_types_T);
  space_T.set_bc_values(bc_values_T);
  space_T.set_uniform_order(P_INIT);
  H1Space space_phi(&mesh_phi, &shapeset);
  space_phi.set_bc_types(bc_types_phi);
  space_phi.set_bc_values(bc_values_phi);
  space_phi.set_uniform_order(P_INIT);

  // assign degrees of freedom
  int ndofs = assign_dofs(2, &space_T, &space_phi);

  // visualisation
  char title[100];
  ScalarView sview_T("Solution computed for the temperature T", 0, 0, 500, 400);
  sview_T.fix_scale_width(80);
  ScalarView sview_phi("Solution computed for the neutron flux phi", 0, 450, 500, 400);
  sview_phi.fix_scale_width(80);
  ScalarView sview_T_exact("Exact solution for the temperature T", 550, 0, 500, 400);
  sview_T_exact.fix_scale_width(80);
  ScalarView sview_phi_exact("Exact solution for the neutron flux phi", 550, 450, 500, 400);
  sview_phi_exact.fix_scale_width(80);
  OrderView ordview_T("Order for the temperature T", 1100, 0, 500, 400);
  ordview_T.fix_scale_width(80);
  OrderView ordview_phi("Order for the neutron flux phi", 1100, 450, 500, 400);
  ordview_phi.fix_scale_width(80);

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof_est, graph_dof_exact, graph_cpu_est, graph_cpu_exact;

  // create a selector which will select optimal candidate
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER, &shapeset);

  // solutions for the Newton's iteration and time stepping
  Solution T_prev_newton, T_prev_time,
           phi_prev_newton, phi_prev_time;

  // initialize the weak formulation
  WeakForm wf(2);
  wf.add_biform(0, 0, jac_TT, jac_TT_ord, H2D_UNSYM, H2D_ANY, 1, &T_prev_newton);
  wf.add_biform(0, 1, jac_Tphi, jac_Tphi_ord, H2D_UNSYM, H2D_ANY, 0);
  wf.add_liform(0, res_T, res_T_ord, H2D_ANY, 3, &T_prev_newton, &T_prev_time, &phi_prev_newton);
  wf.add_biform(1, 0, jac_phiT, jac_phiT_ord, H2D_UNSYM, H2D_ANY, 2, &phi_prev_newton, &T_prev_newton);
  wf.add_biform(1, 1, jac_phiphi, jac_phiphi_ord, H2D_UNSYM, H2D_ANY, 1, &T_prev_newton);
  wf.add_liform(1, res_phi, res_phi_ord, H2D_ANY, 3, &phi_prev_newton, &phi_prev_time, &T_prev_newton);

  // initialize the nonlinear system and solver
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(2, &space_T, &space_phi);
  nls.set_pss(1, &pss);

  // set initial conditions at zero time level
  T_prev_time.set_exact(&mesh_T, T_exact);
  phi_prev_time.set_exact(&mesh_phi, phi_exact);

  // exact solutions for error evaluation
  ExactSolution T_solution(&mesh_T, T_exact),
                phi_solution(&mesh_phi, phi_exact);

  // error
  double T_error,
         phi_error,
         error;

  // time stepping loop
  for (int t_step = 1; TIME<=T_FINAL; t_step++) {
    TIME+=TAU;
    info("-------- Time step %d, t = %g s ------------------------------------------------------", t_step, TIME);
    
    Solution T_coarse, T_fine,
             phi_coarse, phi_fine;

    if (t_step==1) {
      T_fine.copy(&T_prev_time);
      phi_fine.copy(&phi_prev_time);
    }

    // periodic global derefinements
    if (t_step > 1 && t_step % UNREF_FREQ == 0) {
      info("---- Time step %d, global derefinement.", t_step);
      mesh_T.copy(&basemesh);
      mesh_phi.copy(&basemesh);
      space_T.set_uniform_order(P_INIT);
      space_phi.set_uniform_order(P_INIT);
      ndofs = assign_dofs(2, &space_T, &space_phi);

      // project the fine mesh solutions on the globally derefined meshes
      info("---- Time step %d, projecting fine mesh solution on globally derefined mesh.", t_step);
      nls.set_ic(&T_fine, &phi_fine, &T_prev_time, &phi_prev_time, PROJ_TYPE);
    }

    // set initial conditions for the Newton's loop on the coarse mesh
    nls.set_ic(&T_prev_time, &phi_prev_time, &T_prev_newton, &phi_prev_newton, PROJ_TYPE);

    // Newton's loop on the coarse mesh
    info("---- Time step %d, adaptivity step 1, Newton solve on coarse mesh.", t_step);
    if (!nls.solve_newton_2(&T_prev_newton, &phi_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER))
      error("Newton's method did not converge.");

    // store coarse solution
    T_coarse.copy(&T_prev_newton);
    phi_coarse.copy(&phi_prev_newton);

    // adaptivity loop
    double err_est;
    for(int a_step = 1; ; a_step++) {
      // project the fine mesh solution on the new coarse mesh
      if (a_step>1) {
	info("---- Time step %d, adaptivity step %d, projecting fine mesh solution on new coarse mesh.", t_step, a_step);
	nls.set_ic(&T_fine, &phi_fine, &T_coarse, &phi_coarse, PROJ_TYPE);
       }

      info("---- Time step %d, adaptivity step %d, Newton solve on fine mesh.", t_step, a_step);
      // reference nonlinear system
      RefNonlinSystem rnls(&nls);
      rnls.prepare();

      // set initial condition for the Newton's method on the fine mesh
      rnls.set_ic(&T_fine, &phi_fine, &T_prev_newton, &phi_prev_newton, PROJ_TYPE);

      // Newton's loop on the fine mesh
      if (!rnls.solve_newton_2(&T_prev_newton, &phi_prev_newton, NEWTON_TOL_FINE, NEWTON_MAX_ITER))
        error("Newton's method did not converge.");

      // store fine solution
      T_fine.copy(&T_prev_newton);
      phi_fine.copy(&phi_prev_newton);

      // visualization of intermediate solution and mesh during adaptivity
      sprintf(title, "Solution for T, t = %g, adapt step %d", TIME, a_step);
      sview_T.set_title(title);
      sview_T.show(&T_fine);
      sprintf(title, "Solution for phi, t = %g, adapt step %d", TIME, a_step);
      sview_phi.set_title(title);
      sview_phi.show(&phi_fine);
      sprintf(title, "Fine mesh for T, t = %g, adapt step %d", TIME, a_step);
      ordview_T.set_title(title);
      ordview_T.show(rnls.get_space(0));
      sprintf(title, "Fine mesh for phi, t = %g, adapt step %d", TIME, a_step);
      ordview_phi.set_title(title);
      ordview_phi.show(rnls.get_space(1));

      // calculate element errors and total error estimate
      H1Adapt hp(Tuple<Space*>(&space_T, &space_phi));
      hp.set_solutions(Tuple<Solution*>(&T_coarse, &phi_coarse), Tuple<Solution*>(&T_fine, &phi_fine));
      err_est = hp.calc_error() * 100;
      // report error
      info("Error estimate with hp.calc_error_2: %g%%", err_est);

      // compute exact error with function h1_error for comparision with hp.calc_error_2
      T_error = h1_error(&T_coarse, &T_fine) * 100;
      phi_error = h1_error(&phi_coarse, &phi_fine) * 100;
      error = std::max(T_error, phi_error);
      //      info("Error estimate with h1_error for T: %g %%", T_error);
      //      info("Error estimate with h1_error for phi: %g %%", phi_error);
      info("Error estimate with h1_error: %g %%", error);

      // if err_est too large, adapt the mesh
      if (err_est < ERR_STOP) break;
      else {
        hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
        ndofs = assign_dofs(2, &space_T, &space_phi);
        if (ndofs >= NDOF_STOP) break;
      }

      // compute exact solution for comparison with computational results
      T_solution.update(&mesh_T, T_exact);
      phi_solution.update(&mesh_phi, phi_exact);

      // show exact solution
      sprintf(title, "Exact solution for T, t = %g", TIME);
      sview_T_exact.set_title(title);
      sview_T_exact.show(&T_solution);
      sprintf(title, "Exact solution for phi, t = %g", TIME);
      sview_phi_exact.set_title(title);
      sview_phi_exact.show(&phi_solution);

      // moving to new adaptivity step
    }

    // update previous time level solution
    T_prev_time.copy(&T_fine);
    phi_prev_time.copy(&phi_fine);

    // compute exact error
    T_error = h1_error(&T_prev_time, &T_solution) * 100;
    phi_error = h1_error(&phi_prev_time, &phi_solution) * 100;
    error = std::max(T_error, phi_error);
    info("Exact solution error for T (H1 norm): %g %%", T_error);
    info("Exact solution error for phi (H1 norm): %g %%", phi_error);
    info("Exact solution error (maximum): %g %%", error);
  }

  View::wait();

  return 0;
}
