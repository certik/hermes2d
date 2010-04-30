#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "solver_umfpack.h"
#include "function.h"
#include "math.h"
#include <iostream>


// Neutronics/Heat Conduction Test Case
//
// -->PDEs:
//
// 1/v1 d/dt phi_1 = div(D_1 grad phi_1) + (nu Sigma_f1 - Sigma_r1) phi_1 + Sigma_f2 phi_2 
//                        + lambda_1 C1 + lambda_2 C2 + Q_1
//
// 1/v2 d/dt phi_2 = div(D_2 grad phi_2) - Sigma_r2 phi_2 + Sigma_12 phi_1 + Q_2
//
// d/dt C_i = beta_1i nu Sigma_f1 phi_1 + beta_2i nu Sigma_f2 phi_2 - lambda_i C_i    i = 1,2
//
// rho c_p d/dt T = div(k grad T) + kappa_1 Sigma_f1 phi_1 + kappa_2 Sigma_f2 phi_2 + Q_T
//
// -->Domain: rectangle Lx Ly
//
// -->BC: homogeneous Dirichlet
//


double TIME = 0.0;
const int P_INIT = 2;                 // Initial polynomial degree of all mesh elements

const double TAU = 0.1;               // Time step
const double T_FINAL = 10.0;          // Time interval length
const int TIME_MAX_ITER = int(T_FINAL / TAU);

const int PROJ_TYPE = 1;               // For the projection of the initial condition 
                                       // on the initial mesh: 1 = H1 projection, 0 = L2 projection

const int INIT_GLOB_REF_NUM = 4;       // Number of initial uniform mesh refinements
const int INIT_BDY_REF_NUM = 0;        // Number of initial refinements towards boundary

const double NEWTON_TOL = 1e-6;        // Stopping criterion for the Newton's method
const int NEWTON_MAX_ITER = 100;       // Maximum allowed number of Newton iterations

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
const double Tref = 0.0;          // Temperature at boundary

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
  return xsa_ref + doppler_coeff * (sqrt(T) - sqrt(Tref));
}

// Derivative of the removal cross section with respect to temperature
template<typename Real>
Real dxsrem_dT(Real T) {
  return doppler_coeff / (2*sqrt(T+1.0e-10));
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
  //return BC_ESSENTIAL;
  return BC_NATURAL;
}
 

// Dirichlet boundary condition values
scalar bc_values_T(int marker, double x, double y)
{
  return Tref;
}


int bc_types_phi(int marker)
{
  //return BC_ESSENTIAL;
  return BC_NATURAL;
}
 
scalar bc_values_phi(int marker, double x, double y)
{
  return 0.0;
}

/*********************************************************************************************************/
/********** Jacobian and residual  ***********/


// Heat conduction equation
template<typename Real, typename Scalar>
Scalar jac_TT(int n, double *wt, Func<Real> *uj, Func<Real> *ui, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.0;
  Func<Scalar>* T_prev_newton = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (rho * cp * uj->val[i] * ui->val[i] / TAU
		       + dk_dT(T_prev_newton->val[i]) * uj->val[i] * (T_prev_newton->dx[i] * ui->dx[i]
								     + T_prev_newton->dy[i] * ui->dy[i])
                       + k(T_prev_newton->val[i]) * (uj->dx[i] * ui->dx[i]
						     + uj->dy[i] * ui->dy[i]));
  return result;
}

Ord jac_TT_ord(int n, double *wt, Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(10);
}


template<typename Real, typename Scalar>
Scalar jac_Tphi(int n, double *wt, Func<Real> *uj, Func<Real> *ui, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (-kappa * xsfiss * uj->val[i] * ui->val[i]);
  return result;
}

Ord jac_Tphi_ord(int n, double *wt, Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(10);
}

template<typename Real, typename Scalar>
Scalar res_T(int n, double *wt, Func<Real> *ui, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.0;
  Func<Scalar>* T_prev_newton = ext->fn[0];
  Func<Scalar>* T_prev_time = ext->fn[1];
  Func<Scalar>* phi_prev_newton = ext->fn[2];
  for (int i = 0; i < n; i++)
    result += wt[i] * (rho * cp * (T_prev_newton->val[i] - T_prev_time->val[i]) / TAU * ui->val[i]
                       + k(T_prev_newton->val[i]) * (T_prev_newton->dx[i] * ui->dx[i]
						     + T_prev_newton->dy[i] * ui->dy[i])
		       - kappa * xsfiss * phi_prev_newton->val[i] * ui->val[i]
                       - qT(e->x[i], e->y[i]) * ui->val[i]);
  return result;
}

Ord res_T_ord(int n, double *wt, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Neutronics equation
template<typename Real, typename Scalar>
Scalar jac_phiphi(int n, double *wt, Func<Real> *uj, Func<Real> *ui, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.0;
  Func<Scalar>* T_prev_newton = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (invvel * uj->val[i] * ui->val[i] / TAU
                       + xsdiff * (uj->dx[i] * ui->dx[i]
				   + uj->dy[i] * ui->dy[i])
		       + xsrem(T_prev_newton->val[i]) * uj->val[i] * ui->val[i]
		       - nu * xsfiss * uj->val[i] * ui->val[i]);
  return result;
}

Ord jac_phiphi_ord(int n, double *wt, Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(10);
}

template<typename Real, typename Scalar>
Scalar jac_phiT(int n, double *wt, Func<Real> *uj, Func<Real> *ui, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.0;
  Func<Scalar>* phi_prev_newton = ext->fn[0];
  Func<Scalar>* T_prev_newton = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * (dxsrem_dT(T_prev_newton->val[i]) * uj->val[i] * phi_prev_newton->val[i] * ui->val[i]);
  return result;
}

Ord jac_phiT_ord(int n, double *wt, Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(10);
}

template<typename Real, typename Scalar>
Scalar res_phi(int n, double *wt, Func<Real> *ui, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.0;
  Func<Scalar>* phi_prev_newton = ext->fn[0];
  Func<Scalar>* phi_prev_time = ext->fn[1];
  Func<Scalar>* T_prev_newton = ext->fn[2];
  for (int i = 0; i < n; i++)
    result += wt[i] * (invvel * (phi_prev_newton->val[i] - phi_prev_time->val[i]) / TAU * ui->val[i]
                       + xsdiff * (phi_prev_newton->dx[i] * ui->dx[i]
				   + phi_prev_newton->dy[i] * ui->dy[i])
		       + xsrem(T_prev_newton->val[i]) * phi_prev_newton->val[i] * ui->val[i] 
                       - nu * xsfiss * phi_prev_newton->val[i] * ui->val[i]
		       - q(e->x[i], e->y[i]) * ui->val[i]);
  return result;
}

Ord res_phi_ord(int n, double *wt, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}


/*********************************************************************************************************/
/********** Initial conditions ***********/


scalar T_ic(double x, double y, double& dx, double& dy)
{
  dx = 1.0;
  dy = 0.0;
  return x;
}

scalar phi_ic(double x, double y, double& dx, double& dy)
{
  dx = 0.0;
  dy = 1.0;
  return y;
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
  dx = CT*(1+tanh(rT*TIME))*PI/LX*(PI*x/LX)*sin(PI*y/LY);
  dy = CT*(1+tanh(rT*TIME))*sin(PI*x/LX)*PI/LY*sin(PI*y/LY);
  return CT*(1+tanh(rT*TIME))*sin(PI*x/LX)*sin(PI*y/LY);
}

/*********************************************************************************************************/
/********** This is the main ***********/


int main(int argc, char* argv[])
{
  std::cout<<"TIME_MAX_ITER = "<<TIME_MAX_ITER<<std::endl;
  TIME_MAX_ITER;
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // perform initial refinements
  for (int i=0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create H1 spaces
  H1Space space_T(&mesh, &shapeset);
  space_T.set_bc_types(bc_types_T);
  space_T.set_bc_values(bc_values_T);
  space_T.set_uniform_order(P_INIT);
  H1Space space_phi(&mesh, &shapeset);
  space_phi.set_bc_types(bc_types_phi);
  space_phi.set_bc_values(bc_values_phi);
  space_phi.set_uniform_order(P_INIT);

  int ndofs = 0;
  ndofs += space_T.assign_dofs(ndofs);
  ndofs += space_phi.assign_dofs(ndofs);

  // solutions for the Newton's iteration and time stepping
  Solution T_prev_newton, T_prev_time,
           phi_prev_newton, phi_prev_time;

  // exact solutions for error evaluation
  ExactSolution T_solution(&mesh, T_exact),
                phi_solution(&mesh, phi_exact);

  // exact errors
  double T_error,
         phi_error,
         error;

  // setting initial conditions
  T_prev_time.set_exact(&mesh, T_exact);
  //  T_prev_newton.copy(&T_prev_time);
  phi_prev_time.set_exact(&mesh, phi_exact);
  //  phi_prev_newton.copy(&phi_prev_time);

  // visualisation
  ScalarView sview_T("Solution for the temperature T", 0, 0, 500, 400);
  ScalarView sview_phi("Solution for the neutron flux phi", 0, 500, 500, 400);
  ScalarView sview_T_exact("Solution for the temperature T", 550, 0, 500, 400);
  ScalarView sview_phi_exact("Solution for the neutron flux phi", 550, 500, 500, 400);

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
  nls.set_ic(&T_prev_time, &phi_prev_time, &T_prev_newton, &phi_prev_newton, PROJ_TYPE);

  // time stepping loop
  int t_step = 1;
  do {
    TIME += TAU;

    info("!---- Time step %d, t = %g s", t_step, TIME); t_step++;

    // Newton's method
    if (!nls.solve_newton_2(&T_prev_newton, &phi_prev_newton, NEWTON_TOL, NEWTON_MAX_ITER)) error("Newton's method did not converge.");

    // update previous time level solution
    T_prev_time.copy(&T_prev_newton);
    phi_prev_time.copy(&phi_prev_newton);

    // show the new time level solution
    char title[100];
    sprintf(title, "Solution for T, t = %g", TIME);
    sview_T.set_title(title);
    sview_T.show(&T_prev_newton);
    sprintf(title, "Solution for phi, t = %g", TIME);
    sview_phi.set_title(title);
    sview_phi.show(&phi_prev_newton);

    // compute exact solution for comparison with computational results
    // COMMENT: I added membre function update to the class ExactSolution
    //          This is helpful for time-dependent problems
    T_solution.update(&mesh, T_exact);
    phi_solution.update(&mesh, phi_exact);

    // show exact solution
    sprintf(title, "Exact solution for T, t = %g", TIME);
    sview_T_exact.set_title(title);
    sview_T_exact.show(&T_solution);
    sprintf(title, "Exact solution for phi, t = %g", TIME);
    sview_phi_exact.set_title(title);
    sview_phi_exact.show(&phi_solution);

    // compute exact error
    T_error = h1_error(&T_prev_time, &T_solution) * 100;
    phi_error = h1_error(&phi_prev_time, &phi_solution) * 100;
    error = std::max(T_error, phi_error);
    info("Exact solution error for T (H1 norm): %g %%", T_error);
    info("Exact solution error for phi (H1 norm): %g %%", phi_error);
    info("Exact solution error (maximum): %g %%", error);
  }
  while (t_step < TIME_MAX_ITER);

  // wait for all views to be closed
  View::wait();

  return 0;
}
