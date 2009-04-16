#include "hermes2d.h"
#include "solver_umfpack.h"

// 
//  PDE: non-stationary complex Gross-Pitaevski equation
//  describing resonances in Bose-Einstein condensates
//
//  hp-adaptivity with dynamical meshes
//
//  ih \partial \psi/\partial t = -h^2/(2m) \Delta \psi +
//  g \psi |\psi|^2 + 1/2 m \omega^2 (x^2 + y^2) \psi
//
//  Domain: square of edge length 2 centered at origin
//
//  BC:  homogeneous Dirichlet everywhere on th eboundary
//
//  Time-stepping: either implicit Euler or Crank-Nicolson
//

/********** Problem parameters ***********/ 

double H = 1;         //Planck constant 6.626068e-34;
double M = 1; 
double G = 1; 
double OMEGA = 1;     
int TIME_DISCR = 2;    // 1 for implicit Euler, 2 for Crank-Nicolson
int PROJ_TYPE = 1;     // 1 for H1 projections, 2 for L2 projections
double T_FINAL = 2;    // time interval length
double TAU = 0.001;    // time step
int P_INIT = 2;        // initial polynomial degree
int REF_INIT = 3;      // number of initial uniform refinements

/********** Definition of initial conditions ***********/ 

scalar fn_init(double x, double y, scalar& dx, scalar& dy) {
  return complex(exp(-10*(x*x + y*y)), 0);
}

/********** Definition of boundary conditions ***********/ 

int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

complex bc_values(int marker, double x, double y)
{
  return complex(0.0, 0.0);
}

/********** Definition of Jacobian matrices and residual vectors ***********/ 

# include "forms.cpp"

/********** Definition of linear and bilinear forms for Hermes ***********/ 

Solution Psi_prev, // previous time step solution, for the time integration method
         Psi_iter; // solution converging during the Newton's iteration

// Implicit Euler method (1st-order in time)
complex linear_form_0_euler(RealFunction* fv, RefMap* rv)
{ return F_euler(&Psi_prev, &Psi_iter, fv, rv); }
complex bilinear_form_0_0_euler(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J_euler(&Psi_iter, fu, fv, ru, rv); }

// Crank-Nicolson method (2nd-order in time)
complex linear_form_0_cranic(RealFunction* fv, RefMap* rv)
{ return F_cranic(&Psi_prev, &Psi_iter, fv, rv); }
complex bilinear_form_0_0_cranic(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J_cranic(&Psi_iter, fu, fv, ru, rv); }


// *************************************************************

int main(int argc, char* argv[])
{
  Mesh mesh, basemesh;
  basemesh.load("square.mesh");
  for(int i = 0; i < REF_INIT; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);
  
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  space.assign_dofs();    

  WeakForm wf(1);
  if(TIME_DISCR == 1) {
    wf.add_biform(0, 0, bilinear_form_0_0_euler, UNSYM, ANY, 1, &Psi_iter);
    wf.add_liform(0, linear_form_0_euler, ANY, 2, &Psi_iter, &Psi_prev);
  }
  else {
    wf.add_biform(0, 0, bilinear_form_0_0_cranic, UNSYM, ANY, 1, &Psi_iter);
    wf.add_liform(0, linear_form_0_cranic, ANY, 2, &Psi_iter, &Psi_prev);
  }
  
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  char title[100];  
  ScalarView view("", 0, 0, 600, 600);
  OrderView ordview("Polynomial Orders", 0, 700, 600, 600);
  
  GnuplotGraph graph_err;
  graph_err.set_captions("","Time step","Error");
  graph_err.add_row();
  GnuplotGraph graph_dofs;
  graph_dofs.set_captions("","Time step","DOFs");
  graph_dofs.add_row();
  
  // setting initial condition at zero time level
  Psi_prev.set_exact(&mesh, fn_init);
  nls.set_ic(&Psi_prev, &Psi_prev, PROJ_TYPE);
  Psi_iter.copy(&Psi_prev);

  // view initial guess for Newton's method
  // satisfies BC conditions 
  sprintf(title, "Initial iteration");
  view.set_title(title);
  view.show(&Psi_iter);    
  ordview.show(&space);
  //view.wait_for_keypress(); // this may cause graphics problems

  Solution sln, rsln;
  // time stepping
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int n = 1; n <= nstep; n++) 
  {

    info("\n---- Time step %d -----------------------------------------------------------------", n);
    
    // unrefinements
    if (n % 10 == 2) {
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
    }
    
    int at = 0;
    bool done = false;
    double err;
    do 
    {
     info("\n---- Time step %d, adaptivity step %d ---------------------------------------------\n", n, ++at);

    
      int it = 1; 
      double res_l2_norm; 
      //space.assign_dofs();    
      if (n > 1 || at > 1) nls.set_ic(&rsln, &Psi_iter);
      else nls.set_ic(&Psi_iter, &Psi_iter);
      do
      {
        info("\n---- Time step %d, adaptivity step %d, Newton step %d (coarse solution) ----------\n", n, at, it++);
        
        nls.assemble();
        nls.solve(1, &sln);
  
        res_l2_norm = nls.get_residuum_l2_norm(); 
        info("Residuum L2 norm: %g\n", res_l2_norm);
   
        Psi_iter.copy(&sln);          
      }
      while (res_l2_norm > 1e-4);
            
      // reference solution
      it = 1;
      RefNonlinSystem rs(&nls);
      rs.prepare();
      if (n > 1 || at > 1) rs.set_ic(&rsln, &Psi_iter);
      else rs.set_ic(&Psi_iter, &Psi_iter);
      do
      {
        info("\n---- Time step %d, adaptivity step %d, Newton step %d (reference solution) --------\n", n, at, it++);
        
        rs.assemble();
        rs.solve(1, &rsln);

        res_l2_norm = rs.get_residuum_l2_norm(); 
        info("Residuum L2 norm: %g\n", res_l2_norm);
  
        Psi_iter.copy(&rsln);          
      }
      while (res_l2_norm > 1e-4);

      H1OrthoHP hp(1, &space);
      err = hp.calc_error(&sln, &rsln) * 100;
      info("Error: %g", err);
      if (err < 0.001) done = true;
      else hp.adapt(0.3, 1);
      space.assign_dofs();    

      // visualization of intermediate solution 
      // and mesh during adaptivity
      sprintf(title, "Time level %d", n);
      //view.set_min_max_range(-1,1);
      view.set_title(title);
      view.show(&Psi_iter);    
      ordview.show(&space);

    }
    while (!done);
      
    // visualization of solution
    sprintf(title, "Time level %d", n);
    view.set_min_max_range(90,100);
    view.set_title(title);
    view.show(&Psi_iter);    
    ordview.show(&space);
    //view.wait_for_keypress();

    graph_err.add_values(0, n, err);
    graph_err.save("error.txt");
    graph_dofs.add_values(0, n, space.get_num_dofs());
    graph_dofs.save("dofs.txt");
    
    // copying result of the Newton's iteration into Psi_prev
    Psi_prev.copy(&Psi_iter);
  }  

  View::wait();
  return 0;
}
