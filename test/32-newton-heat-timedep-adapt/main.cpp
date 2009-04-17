#include "hermes2d.h"
#include "solver_umfpack.h"

// 
//  PDE: non-stationary heat transfer with nonlinear thermal conductivity
//  HEATCAP*dT/dt - div[lambda(T)grad T] = 0
//
//  hp-adaptivity with dynamical meshes
//
//  Domain: square
//
//  BC:  T = 100 on the left, top and bottom edges
//       dT/dn = 0 on the right edge
//

/********** Problem parameters ***********/ 

int TIME_DISCR = 1;   // 1 for implicit Euler, 2 for Crank-Nicolson
int PROJ_TYPE = 1;    // 1 for H1 projections, 0 for L2 projections
double HEATCAP = 1e6; // heat capacity
double TAU = 0.5;     // time step
double T_FINAL = 60;  // time interval length
int P_INIT = 2;       // initial uniform polynomial order
int REF_INIT = 2;     // number of initial refinements

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

/********** Definition of Jacobian matrices and residual vectors ***********/ 

# include "forms.cpp"

/********** Definition of linear and bilinear forms for Hermes ***********/ 

Solution Tprev, // previous time step solution, for the time integration method
         Titer; // solution converging during the Newton's iteration

// Implicit Euler method (1st-order in time)
scalar linear_form_0_euler(RealFunction* fv, RefMap* rv)
{ return F_euler(&Tprev, &Titer, fv, rv); }
scalar bilinear_form_0_0_euler(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J_euler(&Titer, fu, fv, ru, rv); }

// Crank-Nicolson method (2nd-order in time)
scalar linear_form_0_cranic(RealFunction* fv, RefMap* rv)
{ return F_cranic(&Tprev, &Titer, fv, rv); }
scalar bilinear_form_0_0_cranic(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ return J_cranic(&Titer, fu, fv, ru, rv); }

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
  space.set_uniform_order(2);
  space.assign_dofs();    

  WeakForm wf(1);
  if(TIME_DISCR == 1) {
    wf.add_biform(0, 0, bilinear_form_0_0_euler, UNSYM, ANY, 1, &Titer);
    wf.add_liform(0, linear_form_0_euler, ANY, 2, &Titer, &Tprev);
  }
  else {
    wf.add_biform(0, 0, bilinear_form_0_0_cranic, UNSYM, ANY, 1, &Titer);
    wf.add_liform(0, linear_form_0_cranic, ANY, 2, &Titer, &Tprev);
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
  
  // initial condition at zero time level
  Tprev.set_const(&mesh, 0.0);
  nls.set_ic(&Tprev, &Tprev, PROJ_TYPE);
  Titer.copy(&Tprev);

  // view initial guess for Newton's method
  // satisfies BC conditions 
  sprintf(title, "Initial iteration");
  view.set_title(title);
  view.show(&Titer);    
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
      space.assign_dofs();    
      if (n > 1 || at > 1) nls.set_ic(&rsln, &Titer);
      else nls.set_ic(&Titer, &Titer);
      do
      {
        info("\n---- Time step %d, adaptivity step %d, Newton step %d -----------------------------\n", n, at, it++);
        
        nls.assemble();
        nls.solve(1, &sln);
  
        res_l2_norm = nls.get_residuum_l2_norm(); 
        info("Residuum L2 norm: %g\n", res_l2_norm);
   
        Titer.copy(&sln);          
      }
      while (res_l2_norm > 1e-4);
            
      // reference solution
      it = 1;
      RefNonlinSystem rs(&nls);
      rs.prepare();
      if (n > 1 || at > 1) rs.set_ic(&rsln, &Titer);
      else rs.set_ic(&Titer, &Titer);
      do
      {
        info("\n---- Time step %d, adaptivity step %d, Newton step %d (Reference solution) --------\n", n, at, it++);
        
        rs.assemble();
        rs.solve(1, &rsln);

        res_l2_norm = rs.get_residuum_l2_norm(); 
        info("Residuum L2 norm: %g\n", res_l2_norm);
  
        Titer.copy(&rsln);          
      }
      while (res_l2_norm > 1e-4);

      H1OrthoHP hp(1, &space);
      err = hp.calc_error(&sln, &rsln) * 100;
      info("Error: %g", err);
      if (err < 0.001) done = true;
      else hp.adapt(0.3, 1);

    }
    while (!done);
      
    // visualization of solution on the n-th time level
    sprintf(title, "Time level %d", n);
    view.set_min_max_range(90,100);
    view.set_title(title);
    //view.show(&Titer);    // to see reference solution
    view.show(&sln);        // to see the solution
    ordview.show(&space);   // to see hp-mesh
    //view.wait_for_keypress();

    graph_err.add_values(0, n, err);
    graph_err.save("error.txt");
    graph_dofs.add_values(0, n, space.get_num_dofs());
    graph_dofs.save("dofs.txt");
    
    // copying result of the Newton's iteration into Tprev
    Tprev.copy(&Titer);
  }  

  View::wait();
  return 0;
}
