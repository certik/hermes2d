#include "hermes2d.h"
#include "solver_umfpack.h"

// This example shows how to run automatic hp-adaptivity in Hermes.
// The underlying problem is a simple model of an electrostatic 
// micromotor. We strongly recommend that you run the model 
// with hp-adaptivity first (H_ONLY = 0) and then with h-adaptivity
// (H_ONLY = 1 and P_INIT = 1, H_ONLY = 1 and P_INIT = 2). Compare 
// convergence graphs (produced automatically and easy to visualize 
// using Gnuplot) and CPU times. You'll see a difference!
//
// Additional adaptivity examples on this level are 21-lshape (that
// one has an exact solution) and 22-layer (that one shows multiple-level
// hanging nodes nicely). 
//
// PDE: -div[eps_r grad phi] = 0
//
// BC: phi = 0 V on Gamma_1
//     phi = VOLTAGE on Gamma_2
//     u_2 = 0 on Gamma_1
//     du/dn = 0 on the axis of (planar) symmetry and outer boundary
//     
// The following parameters can be changed: 
// 

int P_INIT = 1;           // initial polynomial degree in mesh
double ERR_STOP = 0.01;   // stopping criterion for hp-adaptivity
                          // (rel. error tolerance between the reference 
                          // and coarse solution in percent)
double THRESHOLD = 0.3;   // error threshold for element refinement
int STRATEGY = 0;         // refinement strategy (0, 1, 2, 3 - see adapt_h1.cpp for explanation)
int H_ONLY = 0;           // if H_ONLY == 0 then full hp-adaptivity takes place, otherwise
                          // h-adaptivity is used. Use this parameter to check that indeed adaptive 
                          // hp-FEM converges much faster than adaptive h-FEM
int NDOF_STOP = 40000;    // adaptivity process stops when the number of degrees of freedom grows over 
                          // this limit. This is mainly to prevent h-adaptivity to go on forever.  
double EPS1 = 1.0;        // permittivity in Omega_1
double EPS2 = 10.0;       // permittivity in Omega_2
double VOLTAGE = 50.0;    // voltage on the stator


scalar bc_values(int marker, double x, double y)
{  
  return (marker == 2) ? VOLTAGE : 0.0;
}

scalar biform1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return EPS1 * int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar biform2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return EPS2 * int_grad_u_grad_v(fu, fv, ru, rv);
}


int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("motor.mesh");

  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  H1Space space(&mesh, &shapeset);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  
  WeakForm wf(1);
  wf.add_biform(0, 0, biform1, SYM, 1);
  wf.add_biform(0, 0, biform2, SYM, 2);
  
  ScalarView sview("Coarse solution", 0, 0, 600, 1000);
  VectorView gview("Gradient", 610, 0, 600, 1000);
  OrderView  oview("Polynomial orders", 1220, 0, 600, 1000);
  gview.set_min_max_range(0.0, 400.0);
  
  Solution sln, rsln;
  UmfpackSolver solver;
  
  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Error Convergence", "Degrees of Freedom", "Error [%]");
  graph.add_row("error estimate", "-", "o");


  int it = 1;
  while (1)
  {
    info("\n---- Iteration %d ---------------------------------------------\n", it++);
    space.assign_dofs();
  
    // solve the coarse problem
    LinSystem ls(&wf, &solver);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln);
    
    // view the solution -- this can be slow; for illustration only
    sview.show(&sln);
    gview.show(&sln, &sln, EPS_NORMAL, FN_DX_0, FN_DY_0);
    oview.show(&space);
    info("");

    // solve the fine (reference) problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &rsln);
    
    // calculate errors and adapt the solution
    H1OrthoHP hp(1, &space);
    double error = hp.calc_error(&sln, &rsln) * 100;

    graph.add_values(0, space.get_num_dofs(), error);
    graph.save("convergence.gp");

    if (error < ERR_STOP || ls.get_num_dofs() >= NDOF_STOP) break;
    hp.adapt(THRESHOLD, STRATEGY, H_ONLY);
  }
  
  // show the fine solution - this is the final result, about
  // one order of magnitude more accurate than the coarse solution
  sview.set_title("Fine solution");
  sview.show(&rsln); 
  gview.show(&rsln, &rsln, EPS_HIGH, FN_DX_0, FN_DY_0);
  
  View::wait();
  return 0;
}
