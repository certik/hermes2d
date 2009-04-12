#include "hermes2d.h"
#include "solver_umfpack.h"

//
//  PDE: non-stationary heat transfer with nonlinear thermal conductivity
//  HEATCAP*dT/dt - div[lambda(T)grad T] = 0
// 
//  Solved using dynamical meshes
//
//  Domain: square
//
//  BC:  T = 100 on the left, top and bottom edges
//       dT/dn = 0 on the right edge
//
//  Time-stepping: either implicit Euler or Crank-Nicolson
//

/********** Problem parameters ***********/ 

int TIME_DISCR = 2;   // 1 for implicit Euler, 2 for Crank-Nicolson
int PROJ_TYPE = 1;    // 1 for H1 projections, 0 for L2 projections
double HEATCAP = 1e6; // heat capacity
double TAU = 0.5;     // time step
int NSTEP = 100;      // number of time steps to do
int SUBDIV_INIT = 2;  // number of levels of initial uniform refinement
int P_INIT = 2;       // initial polynomial degree

// thermal conductivity (temperature-dependent
// for any u, this function has to be  positive in the entire domain!
double lam(double T) { return 1.0 + T*T; } 
double dlam_dT(double T) { return 2*T; }

/********** Definition of boundary conditions ***********/ 

int bc_types(int marker)
{
 if (marker == 4 || marker == 1 || marker == 3) return BC_ESSENTIAL;
//   if (marker == 4) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

scalar bc_values(int marker, double x, double y)
{
 if (marker == 4 || marker == 1 || marker == 3) return 100;
//   if (marker == 4) return -4.0 * sqr(y) + 4.0 * y; 
  else return 0.0;
}

/********** Import of Jacobian matrices and residual vectors ***********/ 

# include "forms.cpp"

/********** Definition of linear and bilinear forms for Hermes *********/ 

Solution Tprev_bas,      // previous time step solution on basic mesh
         Tprev_crs,      // previous time step solution on coarse mesh
         Tprev_ref,      // previous time step solution on reference mesh
         Titer_bas,      // iterated solution on basic mesh
         Titer_crs,      // iterated solution on coarse mesh 
         Titer_ref;      // iterated solution on reference mesh 

// Implicit Euler method (1st-order in time)
scalar linear_form_0_euler(RealFunction* fv, RefMap* rv)
{ return F_euler(&Tprev_bas, &Titer_bas, fv, rv); }
scalar bilinear_form_0_0_euler(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J_euler(&Titer_bas, fu, fv, ru, rv); }

// Crank-Nicolson method (2nd-order in time)
scalar linear_form_0_cranic(RealFunction* fv, RefMap* rv)
{ return F_cranic(&Tprev_bas, &Titer_bas, fv, rv); }
scalar bilinear_form_0_0_cranic(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J_cranic(&Titer_bas, fu, fv, ru, rv); }

/********** Main program ***********************************************/ 

int main(int argc, char* argv[])
{
  // meshes
  Mesh 
    mesh_bas,   // basic 
    mesh_ref,   // reference refined
    mesh_crs;   // reference coarsened
  mesh_bas.load("square.mesh");
  for(int i = 0; i < SUBDIV_INIT; i++) mesh_bas.refine_all_elements();

  // precalculated shapeset
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // registration of weak forms
  WeakForm wf(1);
  if(TIME_DISCR == 1) {
    wf.add_biform(0, 0, bilinear_form_0_0_euler, UNSYM, ANY, 1, &Titer_bas);
    wf.add_liform(0, linear_form_0_euler, ANY, 2, &Titer_bas, &Tprev_bas);
  }
  else {
    wf.add_biform(0, 0, bilinear_form_0_0_cranic, UNSYM, ANY, 1, &Titer_bas);
    wf.add_liform(0, linear_form_0_cranic, ANY, 2, &Titer_bas, &Tprev_bas);
  }
  
  // space for the solution
  H1Space space(&mesh_bas, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  space.assign_dofs();

  // matrix solver
  UmfpackSolver umfpack;

  // nonlinear system
  NonlinSystem nls(&wf, &umfpack);
  nls.set_pss(1, &pss);

  // visualization
  char title[100];  
  ScalarView view("", 0, 0, 600, 600);

  // setting the Dirichlet lift to be the 
  // initial condition on thebasic mesh
  nls.set_spaces(1, &space);  
  Tprev_bas.set_dirichlet_lift(&space, &pss);
  nls.set_ic(&Tprev_bas, &Tprev_bas, PROJ_TYPE);
  Titer_bas.copy(&Tprev_bas);

  // view initial condition
  /*
  sprintf(title, "Initial condition");
  view.set_title(title);
  view.show(&Titer_bas);    
  view.wait_for_keypress();
  */  
  
  // time stepping
  Solution sln, sln_ref, sln_crs;
  for(int n = 1; n <= NSTEP; n++) 
  {

    info("\n---- Time step %d -----------------------------------------------", n);

    // ***** COARSE SOLUTION *****
    // unrefine mesh_bas to mesh_crs, 
    // project solution from mesh_bas to mesh_crs, 
    // solve problem on mesh_crs
    mesh_crs.copy(&mesh_bas);  
    mesh_crs.unrefine_all_elements();
    space.set_mesh(&mesh_crs); 
    space.assign_dofs();
    nls.set_spaces(1, &space);  
    Tprev_crs.copy(&Tprev_bas); 
    nls.set_ic(&Tprev_crs, &Tprev_crs, PROJ_TYPE);
    Titer_crs.copy(&Tprev_crs);
    int it = 1; 
    double res_l2_norm; 
    do
    {
      info("\n---- Solution (coarse): tstep %d, iter %d ---------------------------------\n", n, it++);
      nls.assemble();
      nls.solve(1, &sln);
      res_l2_norm = nls.get_residuum_l2_norm(); 
      info("Residuum L2 norm: %g\n", res_l2_norm);
      Titer_crs = sln;
    }
    while (res_l2_norm > 1e-4);
    Tprev_crs.copy(&Titer_crs);

    // visualization of coarse solution 
    //sprintf(title, "Solution (coarse): Time level %d", n);
    //view.set_title(title);
    //view.show(&Titer_crs);    
    //view.wait_for_keypress();

    // ***** BASIC SOLUTION *****
    // project solution from mesh_crs to mesh_bas, 
    // solve problem on mesh_bas
    space.set_mesh(&mesh_bas); 
    space.assign_dofs();
    nls.set_spaces(1, &space);  
    Tprev_bas.copy(&Tprev_crs); 
    nls.set_ic(&Tprev_bas, &Tprev_bas, PROJ_TYPE);
    it = 1; 
    space.assign_dofs();
    do
    {
      info("\n---- Solution (basic): tstep %d, iter %d ---------------------------------\n", n, it++);
      nls.assemble();
      nls.solve(1, &sln);
      res_l2_norm = nls.get_residuum_l2_norm(); 
      info("Residuum L2 norm: %g\n", res_l2_norm);
      Titer_bas = sln;
    }
    while (res_l2_norm > 1e-4);
    Tprev_bas.copy(&Titer_bas);

    // visualization of solution on the n-th time level
    //sprintf(title, "Solution: Time level %d", n);
    //view.set_title(title);
    //view.show(&Titer);    
    //view.wait_for_keypress();

    // ***** REFERENCE SOLUTION *****
    // refine mesh_bas to mesh_ref, 
    // project solution from mesh_bas to mesh_ref, 
    // solve problem on mesh_ref
    mesh_ref.copy(&mesh_bas);  
    mesh_ref.refine_all_elements();
    space.set_mesh(&mesh_ref); 
    space.assign_dofs();
    nls.set_spaces(1, &space);  
    Tprev_ref.copy(&Tprev_bas); 
    nls.set_ic(&Tprev_ref, &Tprev_ref, PROJ_TYPE);
    Titer_ref.copy(&Tprev_ref);
    it = 1; 
    do
    {
      info("\n---- Solution (fine): tstep %d, iter %d ---------------------------------\n", n, it++);
      nls.assemble();
      nls.solve(1, &sln);
      res_l2_norm = nls.get_residuum_l2_norm(); 
      info("Residuum L2 norm: %g\n", res_l2_norm);
      Titer_crs = sln;
    }
    while (res_l2_norm > 1e-4);
    Tprev_ref.copy(&Titer_ref);

    // visualization of fine solution 
    //sprintf(title, "Solution (fine): Time level %d", n);
    //view.set_title(title);
    //view.show(&Titer_ref);    
    //view.wait_for_keypress();

    

    // TODO - ADAPTIVITY 














  } // end of time stepping loop

  View::wait();
  return 0;
}
