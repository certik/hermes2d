#include "hermes2d.h"
#include "solver_umfpack.h"  // defines the class UmfpackSolver

/*solver settings*/
const int P_INIT = 1;             // initial polynomial degree in mesh
const double THRESHOLD = 0.1;     // error threshold for element refinement
const int STRATEGY = 1;           // refinement strategy (0, 1, 2, 3 - see adapt_h1.cpp for explanation)
const bool H_ONLY = 0;        // if H_ONLY == false then full hp-adaptivity takes place, otherwise
                                  // h-adaptivity is used. Use this parameter to check that indeed adaptive
                                  // hp-FEM converges much faster than adaptive h-FEM
const bool ISO_ONLY = false;      // when ISO_ONLY = true, only isotropic refinements are done,
                                  // otherwise also anisotropic refinements are allowed
const int MESH_REGULARITY = -2;   // specifies maximum allowed level of hanging nodes
                                  // -1 ... arbitrary level hangning nodes
                                  // 1, 2, 3,... means 1-irregular mesh, 2-irregular mesh, etc.
                                  // total regularization (0) is not supported in adaptivity
const double ERR_STOP = 0.01;      // stopping criterion for adaptivity (rel. error tolerance between the
                                  // reference and coarse solution in percent)
const int NDOF_STOP = 50000;      // adaptivity process stops when the number of degrees of freedom grows over
                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.



double mi0=4.0*3.141592654E-7;
double Jext=5000000.0;
double freq=5E3;
double omega=2*3.141592654*freq;


//Fe
int miR2=1000;
double gama2=6E6;
double mi2=mi0*miR2;


int bc_types(int marker)
{
  if (marker==1) {return BC_NATURAL;}
  if (marker==2) {return BC_ESSENTIAL;}
  if (marker==3) {return BC_ESSENTIAL;}
  if (marker==4) {return BC_ESSENTIAL;}
  if (marker==5) {return BC_NONE;}
  if (marker==6) {return BC_NONE;}
}

complex bc_values(int marker, double x, double y)
{
  return complex(0.0,0.0);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_Gamma_1(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return 0.0;
}


template<typename Real, typename Scalar>
Scalar bilinear_form_Fe(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  scalar ii = complex(0.0, 1.0);
  return 1/mi2*int_grad_u_grad_v<Real, Scalar>(n, wt, u, v)-ii*omega*gama2*int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_vodic(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (1/mi0*int_grad_u_grad_v<Real, Scalar>(n, wt, u, v));
}

template<typename Real, typename Scalar>
Scalar linear_form_vodic(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return Jext*int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_okoli(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return 1/mi0*(int_grad_u_grad_v<Real, Scalar>(n, wt, u, v));
}


Solution slnF,rslnF;

void elMagPole()
{
  int it = 1;
  bool done = false;

  // load the mesh file
  Mesh mesh;
  mesh.load("domain2.mesh");

  // initial mesh refinement (here you can apply arbitrary
  // other initial refinements, see example 01)

//  mesh.refine_all_elements();          // refines all elements
//  mesh.refine_all_elements();          // refines all elements
// mesh.refine_all_elements();          // refines all elements
// mesh.refine_all_elements();          // refines all elements
//  mesh.refine_all_elements();          // refines all elements

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);
  
  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);


  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0,callback(bilinear_form_Fe),UNSYM,3);
  wf.add_biform(0, 0,callback(bilinear_form_vodic),UNSYM,2);
  wf.add_biform(0, 0,callback(bilinear_form_okoli),UNSYM,1);
  wf.add_liform(0, callback(linear_form_vodic),2);
  wf.add_liform_surf(0, callback(linear_form_surf_Gamma_1), 1);


  ScalarView view("Vector potential A ");
  OrderView  oview("Polynomial orders", 1220, 0, 600, 1000);

  do
  {
  info("\n---- Iteration %d ---------------------------------------------\n", it++);
  space.assign_dofs();
  // initialize the linear system and solver
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);
  
  // assemble the stiffness matrix and solve the system

  sys.assemble();
  sys.solve(1, &slnF);
  // visualize the solution
  view.show(&slnF,EPS_LOW);
  oview.show(&space);

    RefSystem rs(&sys);
    rs.assemble();
    rs.solve(1, &rslnF);

    // calculate element errors and total error estimate
    H1OrthoHP hp(1, &space);
    double err_est = hp.calc_error(&slnF, &rslnF) * 100;
    info("Error estimate: %g%%", err_est);

    // adaptivity step
    if (err_est < ERR_STOP || sys.get_num_dofs() >= NDOF_STOP) done = true;
    else hp.adapt(THRESHOLD, STRATEGY, H_ONLY, ISO_ONLY, MESH_REGULARITY);

  }
  while (done == false);

  View::wait();
}

int main(int argc, char* argv[])
{
  elMagPole();
return 0;
}