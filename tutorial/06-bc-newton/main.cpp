#include "hermes2d.h"
#include "solver_umfpack.h"

// This example is a continuation of examples 03, 04 and 05. It explains
// how to use Newton boundary conditions, and again it shows how a Filter
// is used to visualize the solution gradient (which is infinite at the
// re-entrant corner again)
//
// PDE: Laplace equation -Laplace u = 0 (this equation describes, among
// many other things, also stationary heat transfer in a homogeneous linear
// material).
//
// BC: u = T1 ... fixed temperature on Gamma_3 (Dirichlet)
//     du/dn = 0 ... insulated wall on Gamma_2 and Gamma_4 (Neumann)
//     du/dn = H*(u - T0) ... heat flux on Gamma_1 (Newton)
//
// (Note that the last BC can be written in the form  du/dn - H*u = -h*T0 )
//
// You can play with the parameters below:
//

double T1 = 30.0;            // prescribed temperature on Gamma_3
double T0 = 20.0;            // outer temperature on Gamma_1
double H  = 0.05;            // heat flux on Gamma_1
int P_INIT = 6;              // uniform polynomial degree in the mesh
int UNIFORM_REF_LEVEL = 2;   // number of initial uniform mesh refinements
int CORNER_REF_LEVEL = 12;   // number of mesh refinements towards the re-entrant corner


int bc_types(int marker)
  { return (marker == 3) ? BC_ESSENTIAL : BC_NATURAL; }

scalar bc_values(int marker, double x, double y)
  { return (marker == 3) ? T1 : 0.0; }


scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv); }

scalar bilinear_form_surf_Gamma_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv, EdgePos* ep)
  { return H * surf_int_u_v(fu, fv, ru, rv, ep); }

scalar linear_form_surf_Gamma_1(RealFunction* fv, RefMap* rv, EdgePos* ep)
  { return T0 * H * surf_int_v(fv, rv, ep); }


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");
  for(int i=0; i<UNIFORM_REF_LEVEL; i++) mesh.refine_all_elements();
  mesh.refine_towards_vertex(3, CORNER_REF_LEVEL);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  space.assign_dofs();

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form);
  wf.add_biform_surf(0, 0, bilinear_form_surf_Gamma_1, 1);
  wf.add_liform_surf(0, linear_form_surf_Gamma_1, 1);

  // initialize the linear system and solver
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);

  // assemble the stiffness matrix and solve the system
  Solution sln;
  sys.assemble();
  sys.solve(1, &sln);

  // visualize the solution
  ScalarView view("Solution");
  view.show(&sln);

  // compute and show gradient magnitude
  // (note that the infinite gradient at the re-entrant
  // corner will be truncated for visualization purposes)
  ScalarView gradview("Gradient");
  MagFilter grad(&sln, &sln, FN_DX, FN_DY);
  gradview.show(&grad);

  // waiting for keyboard or mouse input
  View::wait();
  return 0;
}
