#include "hermes2d.h"
#include "solver_umfpack.h"
#include "function.h"

// This class performs the Newton's method 
// Y_{n+1} = Y_n - NEWTON_ALPHA*J^{-1}(Y_n)*F(Y_n)
// for the vector-valued system F(Y) = 0.
// Here, F(Y) is the nonlinear function and 
// J(Y) = DF/DY its Jacobi matrix. The function
// F(Y) is defined analogously to the right-hand
// side in the LinSystem class, and the matrix
// J(Y) analogously to the stiffness matrix.

double RESIDUUM_NORM_TOL = 1e-5;  // The Newton's iteration stops if the norm of the 
                                  // vector F(Y_{n+1}) is less than this value. 
                                  // Currently, the maximum norm is used.
double NEWTON_ALPHA = 1.0;        // Parameter for the Newton's method (0 < NEWTON_ALPHA <= 1)
                                  // Default value is NEWTON_ALPHA = 1.
double TEMP_BDY = 100.0;          // This temperature corresponds to the bdy marker 5,
                                  // other bdy markers have zero Neumann condition

// thermal conductivity (temperature-dependent)
// NOTE: for any T, this function has to be  
// positive in the entire domain!
double lam(double T) { return 1 + pow(T,2); } 
double dlam_dT(double T) { return 2*pow(T,1); }
                                                                                                                      
// definition of boundary conditions
int T_bc_type(int marker)
  { return (marker == 5) ? BC_ESSENTIAL : BC_NATURAL; }
double T_bc_value(int marker, double x, double y)
  { return TEMP_BDY; }

// solution from the previous iteration
Solution Tprev;

inline double J(RealFunction* Tprev, RealFunction* fu, 
                RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = Tprev->get_fn_order() + fu->get_fn_order() 
          + fv->get_fn_order() + ru->get_inv_ref_order() + 4;
  limit_order(o);
  Tprev->set_quad_order(o, FN_DEFAULT);
  fu->set_quad_order(o, FN_DEFAULT);
  fv->set_quad_order(o, FN_DEFAULT);

  double* Tprev_val = Tprev->get_fn_values();
  double* uval = fu->get_fn_values();
  double* vval = fv->get_fn_values();
  double *dTprev_dx, *dTprev_dy, *dudx, *dudy, *dvdx, *dvdy;
  Tprev->get_dx_dy_values(dTprev_dx, dTprev_dy);
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);

  double result0 = 0.0; 
  double3* pt = quad->get_points(o); 
  int np = quad->get_num_points(o); 
  double2x2 *mv, *mu; 
  mu = ru->get_inv_ref_map(o); 
  mv = rv->get_inv_ref_map(o); 
  double* jac = ru->get_jacobian(o); 
  for (int i = 0; i < np; i++, mu++, mv++) {
    //    result0 += pt[i][2] * jac[i] * (
    //  dlam_dT(Tprev_val[i]) * vval[i] * (dTprev_dx[i]*t_dudx + dTprev_dy[i]*t_dudy)
    //  + lam(Tprev_val[i]) * (t_dvdx*t_dudx + t_dvdy*t_dudy)); 
    result0 += pt[i][2] * jac[i] * (
      dlam_dT(Tprev_val[i]) * uval[i] * (dTprev_dx[i]*t_dvdx + dTprev_dy[i]*t_dvdy)
      + lam(Tprev_val[i]) * (t_dudx*t_dvdx + t_dudy*t_dvdy)); 
  }

  return result0;
}

inline double F(RealFunction* Tprev, RealFunction* fu, RefMap* ru)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = Tprev->get_fn_order() + fu->get_fn_order() 
          + ru->get_inv_ref_order() + 4;
  limit_order(o);
  Tprev->set_quad_order(o, FN_DEFAULT);
  fu->set_quad_order(o, FN_DEFAULT);

  double* Tprev_val = Tprev->get_fn_values();
  double* vval = fu->get_fn_values();
  double *dTprev_dx, *dTprev_dy, *dudx, *dudy;
  Tprev->get_dx_dy_values(dTprev_dx, dTprev_dy);
  fu->get_dx_dy_values(dudx, dudy);

  double result = 0.0; 
  double3* pt = quad->get_points(o); 
  int np = quad->get_num_points(o); 
  double2x2 *mv, *mu; 
  mu = ru->get_inv_ref_map(o); 
  // mv = rv->get_inv_ref_map(o); 
  double* jac = ru->get_jacobian(o); 
  for (int i = 0; i < np; i++, mu++, mv++) {
    result += pt[i][2] * jac[i] * (
      lam(Tprev_val[i]) * (dTprev_dx[i]*t_dudx + dTprev_dy[i]*t_dudy)); 
  }
  return result;
}

scalar bilinear_form_sym_0_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J(&Tprev, fu, fv, ru, rv); }

scalar linear_form_0(RealFunction* fv, RefMap* rv)
{ return F(&Tprev, fv, rv); }

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  // mesh.load("hole_new.mesh");
  mesh.load("square.mesh");
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();

  // mesh.refine_towards_boundary(5, 0);

  H1ShapesetBeuchler shapeset;
  PrecalcShapeset pss(&shapeset);

  // space for the temperature
  H1Space T(&mesh, &shapeset);

  // initialize boundary conditions
  T.set_bc_types(T_bc_type);
  T.set_bc_values(T_bc_value);

  // set element polynomial degrees
  T.set_uniform_order(1);

  // assign degrees of freedom and calculate problem size
  int ndofs = 0;
  ndofs += T.assign_dofs(ndofs);

  // setting zero initial condition
  // Tprev.set_zero(&mesh);
  scalar *y_temp = new scalar[ndofs];
  for (int i = 0; i<ndofs; i++) y_temp[i] = 0;
  Tprev.set_fe_solution(&T, &pss, y_temp);
  delete [] y_temp;

  // visualization
  ScalarView Tview("Temperature", 0, 500, 1200, 470);
  // Tview.set_min_max_range(0, TEMP_BDY);
  // Fview.show_scale(false);

  // set up weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form_sym_0_0, SYM, 0, 1, &Tprev);
  wf.add_liform(0, linear_form_0, 0, 1, &Tprev);

  // set up the nonlinear system
  UmfpackSolver umfpack;
  NonlinSystem nsys(&wf, &umfpack);
  nsys.set_spaces(1, &T);
  nsys.set_pss(1, &pss);
  nsys.set_alpha(NEWTON_ALPHA);

  // defining zero initial condition for the Newton's 
  // iteration (fixme - Tprev.set_zero(&mesh) used above 
  // should be enough to do this.
  nsys.set_vec_zero();

  // main loop
  int count = 0;
  double residuum_norm = 1e10;
  do 
  {
    Solution Tsln;
    info("\n*** Iteration %d ***", ++count);
    
    // assemble
    nsys.assemble();

    /* DEBUG
    // output matrix and RHS in Matlab format
    char *filename = new char[20];
    strcpy(filename, "matrix.txt");
    char *varname = new char[20];
    strcpy(varname, "matrix");
    nsys.save_matrix_matlab(filename, varname);
    strcpy(filename, "rhs.txt");
    strcpy(varname, "rhs");
    nsys.save_rhs_matlab(filename, varname);
    */
 
    // solve
    nsys.solve(1, &Tsln);

    // calculate the residuum
    residuum_norm = nsys.get_residuum_max_norm();
    info("Residuum max norm = %g\n", residuum_norm);
    
    // visualization of intermediate Newton's iterations
    //Tview.show(&Tsln, 2*EPS_LOW);    
    //Tview.wait();

    // uncomment one of the following lines to generate a series of video frames
    // vview.save_numbered_screenshot("velocity%03d.bmp", i, true);
    // pview.save_numbered_screenshot("pressure%03d.bmp", i, true);
    // the frames can then be converted to a video file with the command
    // mencoder "mf://velocity*.bmp" -mf fps=20 -o velocity.avi -ovc lavc -lavcopts vcodec=mpeg4
    
    Tprev = Tsln;

    info("Residuum norm = %g\n", residuum_norm);
  } while(residuum_norm > RESIDUUM_NORM_TOL);

  View::wait();
}
