#include "hermes2d.h"
#include "solver_umfpack.h"
#include "_hermes2d_api.h"
#include "forms.h"

// The following parameters can be played with:

const double FINAL_TIME = 3600*72;    // length of time interval
const int P_INIT_w0 = 1;       // polynomial degree for pressure
const int P_INIT_VEL = 1;            // polynomial degree for velocity components
const int P_INIT_w4 = 1;       // polynomial degree for pressure

// global time variable
double TIME = 0;

static void calc_pressure_func(int n,
        scalar* w0,
        scalar* w1,
        scalar* w3,
        scalar* w4,
        scalar* result)
{
  for (int i = 0; i < n; i++)
    result[i] = R/c_v * (w4[i] - (w1[i]*w1[i] + w3[i]*w3[i])/(2*w0[i]));
}

class CalcPressure : public SimpleFilter
{
public:
    CalcPressure(
        MeshFunction* sln1,
        MeshFunction* sln2,
        MeshFunction* sln3,
        MeshFunction* sln4,
        int item1 = FN_VAL,
        int item2 = FN_VAL,
        int item3 = FN_VAL,
        int item4 = FN_VAL
        )
        : SimpleFilter(calc_pressure_func, sln1, sln2, sln3, sln4,
                item1, item2, item3, item4)
    {}
};



int main(int argc, char* argv[])
{
  // import hermes2d

  printf("Importing hermes1d\n");
  // Initialize Python
  Py_Initialize();
  PySys_SetArgv(argc, argv);
  if (import_hermes2d___hermes2d())
      throw std::runtime_error("hermes2d failed to import.");
  cmd("print 'Python initialized'");

  // load the mesh file
  Mesh mesh;
  mesh.load("domain-quad.mesh"); // unstructured triangular mesh available in domain-tri.mesh

  // a-priori mesh refinements
  //mesh.refine_all_elements();
  //mesh.refine_all_elements();
  //mesh.refine_all_elements();
  mesh.refine_towards_boundary(marker_bottom, 4);
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements(2);
  mesh.refine_all_elements(2);
  //mesh.refine_all_elements();
  //mesh.refine_all_elements();
  //mesh.refine_towards_boundary(1, 4);
  //mesh.refine_towards_boundary(3, 4);

  // display the mesh
  //MeshView mview("Navier-Stokes Example - Mesh", 100, 100, 1100, 400);
  //mview.show(&mesh);
  //mview.wait_for_keypress();

  // initialize the shapesets and the cache
  H1Shapeset shapeset_h1;
  PrecalcShapeset pss_h1(&shapeset_h1);

  // H1 spaces for velocities and L2 for pressure
  H1Space s0(&mesh, &shapeset_h1);
  H1Space s1(&mesh, &shapeset_h1);
  H1Space s3(&mesh, &shapeset_h1);
  H1Space s4(&mesh, &shapeset_h1);

  register_bc(s0, s1, s3, s4);

  // set velocity and pressure polynomial degrees
  s0.set_uniform_order(P_INIT_w0);
  s1.set_uniform_order(P_INIT_VEL);
  s3.set_uniform_order(P_INIT_VEL);
  s4.set_uniform_order(P_INIT_w4);

  // assign degrees of freedom
  int ndofs = 0;
  ndofs += s0.assign_dofs(ndofs);
  ndofs += s1.assign_dofs(ndofs);
  ndofs += s3.assign_dofs(ndofs);
  ndofs += s4.assign_dofs(ndofs);

  // initial condition: xprev and yprev are zero
  Solution w0_prev, w1_prev, w3_prev, w4_prev;
  set_ic(mesh, w0_prev, w1_prev, w3_prev, w4_prev);

  // set up weak formulation
  WeakForm wf(4);
  register_forms(wf, w0_prev, w1_prev, w3_prev, w4_prev);

  // visualization
  VectorView w13_view("Current Density [m/s]", 0, 0, 1500, 470);
  ScalarView w0_view("Mass Density [Pa]", 0, 530, 1500, 470);
  ScalarView w4_view("Energy [Pa]", 0, 530, 1500, 470);
  //w13_view.set_min_max_range(0, 2);
  w0_view.show_mesh(false);
  w4_view.show_mesh(false);
  // fixing scale width (for nicer videos). Note: creation of videos is
  // discussed in a separate example
  //vview.fix_scale_width(5);
  //pview.fix_scale_width(5);

  // set up the linear system
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(4, &s0, &s1, &s3, &s4);
  sys.set_pss(4, &pss_h1, &pss_h1, &pss_h1, &pss_h1);
  //sys.set_pss(1, &pss_h1);


  // main loop
  char title[100];
  int num_time_steps = FINAL_TIME / TAU;
  for (int i = 1; i <= num_time_steps; i++)
  {
    TIME += TAU;
    set_iteration(i);

    cmd("print 'Iteration'");
    info("\n---- Time step %d, time = %g -----------------------------------", i, TIME);

    // this is needed to update the time-dependent boundary conditions
    ndofs = 0;
    ndofs += s0.assign_dofs(ndofs);
    ndofs += s1.assign_dofs(ndofs);
    ndofs += s3.assign_dofs(ndofs);
    ndofs += s4.assign_dofs(ndofs);

    // assemble and solve
    Solution w0_sln, w1_sln, w3_sln, w4_sln;
    sys.assemble();
    sys.solve(4, &w0_sln, &w1_sln, &w3_sln, &w4_sln);

    // visualization
    sprintf(title, "Current density, time %g", TIME);
    w13_view.set_title(title);
    w13_view.show(&w1_prev, &w3_prev, EPS_LOW);
    sprintf(title, "Mass Density, time %g", TIME);
    w0_view.set_title(title);
    w0_view.show(&w0_sln);
    sprintf(title, "Pressure, time %g", TIME);
    CalcPressure pressure(&w0_sln, &w1_sln, &w3_sln, &w4_sln);

    w4_view.set_title(title);
    w4_view.show(&pressure);

    w0_prev = w0_sln;
    w1_prev = w1_sln;
    w3_prev = w3_sln;
    w4_prev = w4_sln;
  }

  View::wait();
}
