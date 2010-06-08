#include "hermes2d.h"
#include "solver_umfpack.h"
#define H2D_REPORT_FILE "test.log"

/*
 * This example is only for showing how to use the values from neighbor in linear surface forms.
 * Solution of the original problem is function x^3 + y^3. As simulated input from previous step (the need to have values from neighbor is originally from time dependent
 * problem) is taken exact solution.
 * The example doesn't have any real meaning.
 */


const int P_INIT = 3;          // Polynomial degree of mesh elements
const int INIT_REF_NUM = 3;    // Number of initial uniform mesh refinements

scalar proj_func(double x, double y, double &dx, double &dy)
{
  return x*x*x + y*y*y;
}

// bilinear and linear form defining the projection
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_u_v<Real, Scalar>(n, wt, u, v);
}

// return the value \int v dx
template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++){
    result += wt[i] * ((pow(e->x[i], 1) + pow(e->y[i], 1)) * v->val[i]);
  }
  return result;
}


template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  // here to get values from neighbor you have to use method get_fn_neighbor(). This is for safety.
  for (int i = 0; i < n; i++){
     result += 0.5*wt[i] * (ext->fn[0]->val[i] + ext->get_fn_neighbor(0)->val[i]) * v->val[i];
  }
  return result;
}

// boundary conditions
int bc_types(int marker)
{
   return BC_NONE;
}

int main(int argc, char* argv[])
{
  if (argc < 2) error("Missing mesh file name parameter.");

  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load(argv[1], &mesh);

  // uniform mesh refinements
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  mesh.refine_element(29);
  mesh.refine_element(87);

  // display the mesh
	 MeshView mview("info_neighbor", 100, 100, 500, 500);
	 mview.show(&mesh);
  // wait for keyboard or mouse input
   View::wait("Waiting for keyboard or mouse input.");


  // initialize the shapeset and the cache
  L2Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  // create the L2 space
  L2Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);

  // set uniform polynomial degrees
  space.set_uniform_order(P_INIT);

  // enumerate basis functions
  int ndof = assign_dofs(&space);

 /// BaseView bview;
//  bview.show(&space);
 // bview.wait_for_close();

  Solution sln;
  Solution xprev;
  xprev.set_exact(&mesh, &proj_func);
  // matrix solver
  UmfpackSolver umfpack;

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form));
  wf.add_liform(0, callback(linear_form));
  // if you want to work in linear surface form with values from neighbors use flag H2D_ANY_INNER_EDGE.
  wf.add_liform_surf(0, callback(linear_form_surf), H2D_ANY_INNER_EDGE, 1, &xprev);

  // assemble and solve the finite element problem
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);
  sys.assemble();
  sys.solve(1, &sln);

  // visualize the solution
  ScalarView view1("Solution");
  view1.show(&sln);

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}

