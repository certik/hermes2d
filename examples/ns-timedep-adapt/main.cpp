#include "hermes2d.h"
#include "solver_umfpack.h"

// The time-dependent laminar incompressible Navier-Stokes equations are
// discretized in time via the implicit Euler method. The convective term
// is linearized simply by replacing the velocity in front of the nabla
// operator with the velocity from last time step.
//
// PDE: incompressible Navier-Stokes equations in the form
// \partial v / \partial t - \Delta v / Re + (v \cdot \nabla) v + \nabla p = 0,
// div v = 0
//
// BC: u_1 is a time-dependent constant and u_2 = 0 on Gamma_4 (inlet)
//     u_1 = u_2 = 0 on Gamma_1 (bottom), Gamma_3 (top) and Gamma_5 (obstacle)
//     "do nothing" on Gamma_2 (outlet)
//
// TODO: Implement Crank-Nicolson so that comparisons with implicit Euler can be made
//
// The following parameters can be changed:
//

const double RE = 200.0;             // Reynolds number
const double VEL_INLET = 1.0;        // inlet velocity (reached after STARTUP_TIME)
const double STARTUP_TIME = 1.0;     // during this time, inlet velocity increases gradually
                                     // from 0 to VEL_INLET, then it stays constant
const double TAU = 0.1;              // time step
const double FINAL_TIME = 20.0;      // length of time interval
const int P_INIT_VEL = 2;            // initial polynomial degree for velocity components
const int P_INIT_PRESSURE = 1;       // initial polynomial degree for pressure
                                     // Note: P_INIT_VEL should always be greater than
                                     // P_INIT_PRESSURE because of the inf-sup condition
const double H = 5;                  // domain height (necessary to define the parabolic
                                     // velocity profile at inlet)

//  to better understand boundary conditions
const int marker_bottom = 1;
const int marker_right  = 2;
const int marker_top = 3;
const int marker_left = 4;
const int marker_obstacle = 5;

// global time variable
double current_time = 0;

// adaptivity
const double space_tol = 0.5;
const double thr = 0.3;
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// definition of boundary conditions
int xvel_bc_type(int marker) {
  if (marker == marker_right) return BC_NONE;
  else return BC_ESSENTIAL;
}

int yvel_bc_type(int marker) {
  if (marker == marker_right) return BC_NONE;
  else return BC_ESSENTIAL;
}

int press_bc_type(int marker)
  { return BC_NONE; }

scalar xvel_bc_value(int marker, double x, double y) {
  if (marker == marker_left) {
    // time-dependent inlet velocity
    //double val_y = VEL_INLET; //constant profile
    double val_y = VEL_INLET * y*(H-y) / (H/2.)/(H/2.); //parabolic profile with peak VEL_INLET at y = H/2
    if (current_time <= STARTUP_TIME) return val_y * current_time/STARTUP_TIME;
    else return val_y;
  }
  else return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar bilinear_form_sym_0_0_1_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) / RE + int_u_v<Real, Scalar>(n, wt, u, v) / TAU;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_unsym_0_0_1_1(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_w_nabla_u_v<Real, Scalar>(n, wt, ext->fn[0], ext->fn[1], u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_u_v<Real, Scalar>(n, wt, ext->fn[0], v) / TAU;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_unsym_0_2(int n, double *wt, Func<Real> *p, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return - int_u_dvdx<Real, Scalar>(n, wt, p, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_unsym_1_2(int n, double *wt, Func<Real> *p, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return - int_u_dvdy<Real, Scalar>(n, wt, p, v);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Magnitude of velocity filter

void mag(int n, scalar* a, scalar* dadx, scalar* dady,
         scalar* b, scalar* dbdx, scalar* dbdy,
         scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = sqrt(sqr(a[i]) + sqr(b[i]));
    outdx[i] = (0.5 / out[i]) * (2.0 * a[i] * dadx[i] + 2.0 * b[i] * dbdx[i]);
    outdx[i] = (0.5 / out[i]) * (2.0 * a[i] * dady[i] + 2.0 * b[i] * dbdy[i]);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////   ADAPTIVITY   ////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

static double* cmp_err;
static int compare(const void* p1, const void* p2)
{
  const int (*e1) = ((const int*) p1);
  const int (*e2) = ((const int*) p2);
  return cmp_err[(*e1)] < cmp_err[(*e2)] ? 1 : -1;
}

// calculates element errors between 2 solutions (coarse and normal, normal and fine)
double calc_error(MeshFunction* sln, MeshFunction* rsln, int*& esort, double*& errors)
{
  int i, j;

  Mesh* cmesh = sln->get_mesh();
  int max = cmesh->get_max_element_id();
  errors = new double[max];
  memset(errors, 0, sizeof(double) * max);
  int nact = cmesh->get_num_active_elements();
  esort = new int[nact];

  double total_error = 0.0, total_norm = 0.0;

  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  rsln->set_quad_2d(quad);

  Mesh* meshes[2] = { sln->get_mesh(), rsln->get_mesh() };
  Transformable* tr[2] = { sln, rsln };
  Traverse trav;
  trav.begin(2, meshes, tr);

  Element** ee;
  while ((ee = trav.get_next_state(NULL, NULL)) != NULL)
  {
    update_limit_table(ee[0]->get_mode());

    RefMap* crm = sln->get_refmap();
    RefMap* frm = rsln->get_refmap();

    double err  = int_l2_error(sln, rsln, crm, frm);
    total_norm += int_l2_norm (rsln, frm);

    errors[ee[0]->id] += err;
    total_error += err;
  }
  trav.finish();

  Element* e;
  for_all_inactive_elements(e, cmesh)
    errors[e->id] = -1.0;

  int k = 0;
  for_all_active_elements(e, cmesh)
  {
    errors[e->id] /= total_norm;
    esort[k++] = e->id;
  }

  cmp_err = errors;
  qsort(esort, nact, sizeof(int), compare);

  return sqrt(total_error / total_norm);
}


// refines or coarses the mesh using element errors calculated in calc_error()
void adapt_mesh(bool& done, Mesh* mesh, Mesh* cmesh, Space* space, int* esort0, double* errors0, int* esort, double* errors)
{
  int i, j;
  double err0 = esort[0];
  if (!done) // refinements
  {
    for (i = 0; i < mesh->get_num_active_elements(); i++)
    {
      int id = esort[i];
      double err = errors[id];

      if (err < thr * errors[esort[0]]) { break; }
      err0 = err;

      Element* e;
      e = mesh->get_element(id);
      mesh->refine_element(id);
      for (j = 0; j < 4; j++)
        space->set_element_order(e->sons[j]->id, space->get_element_order(id));
    }
  }

  if (done)  // coarsening
  {
    for (i = 0; i < cmesh->get_num_active_elements(); i++)
    {
      int id = esort0[i];
      double err = errors0[id];

      if (err < thr * errors[esort[0]])
      {
        Element* e;
        e = mesh->get_element(id);
        if (!(e->active))
        {
          int o = space->get_element_order(e->sons[0]->id);
          mesh->unrefine_element(id);
          space->set_element_order(id, o);
        }
      }
    }
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);
  mesh.refine_towards_boundary(1, 2);
  mesh.refine_towards_boundary(3, 2);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // spaces for velocities and pressure
  H1Space xvel(&mesh, &shapeset);
  H1Space yvel(&mesh, &shapeset);
  H1Space press(&mesh, &shapeset);

  // initialize boundary conditions
  xvel.set_bc_types(xvel_bc_type);
  xvel.set_bc_values(xvel_bc_value);
  yvel.set_bc_types(yvel_bc_type);
  press.set_bc_types(press_bc_type);

  // set velocity and pressure polynomial degrees
  xvel.set_uniform_order(P_INIT_VEL);
  yvel.set_uniform_order(P_INIT_VEL);
  press.set_uniform_order(P_INIT_PRESSURE);

  // assign degrees of freedom
  int ndofs = 0;
  ndofs += xvel.assign_dofs(ndofs);
  ndofs += yvel.assign_dofs(ndofs);
  ndofs += press.assign_dofs(ndofs);

  // initial BC: xprev and yprev are zero
  Solution xprev, yprev;
  xprev.set_zero(&mesh);
  yprev.set_zero(&mesh);

  // set up weak formulation
  WeakForm wf(3);
  wf.add_biform(0, 0, callback(bilinear_form_sym_0_0_1_1), SYM);
  wf.add_biform(0, 0, callback(bilinear_form_unsym_0_0_1_1), UNSYM, ANY, 2, &xprev, &yprev);
  wf.add_biform(1, 1, callback(bilinear_form_sym_0_0_1_1), SYM);
  wf.add_biform(1, 1, callback(bilinear_form_unsym_0_0_1_1), UNSYM, ANY, 2, &xprev, &yprev);
  wf.add_biform(0, 2, callback(bilinear_form_unsym_0_2), ANTISYM);
  wf.add_biform(1, 2, callback(bilinear_form_unsym_1_2), ANTISYM);
  wf.add_liform(0, callback(linear_form), ANY, 1, &xprev);
  wf.add_liform(1, callback(linear_form), ANY, 1, &yprev);

  // visualization
  VectorView vview("velocity [m/s]", 0, 0, 1500, 470);
  ScalarView pview("pressure [Pa]", 0, 530, 1500, 470);
  vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 0.9);
  pview.show_mesh(false);
  pview.fix_scale_width(80);

  // set up the linear system
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(3, &xvel, &yvel, &press);
  sys.set_pss(1, &pss);

  // main loop
  char title[100];
  for (int i = 1; current_time < FINAL_TIME; i++)
  {


    info("\n---- Time step %d, time = %g -----------------------------------", i, current_time += TAU);

    Solution xcrs, ycrs, pcrs;
    Solution xsln, ysln, psln;
    Solution xref, yref, pref;

    bool done = false;
    int at = 0;
    do
    {
      info("\n*** Adaptive iteration %d ***\n", at++);

      // solve problem
      int ndofs = 0;
      ndofs += xvel.assign_dofs(ndofs);
      ndofs += yvel.assign_dofs(ndofs);
      ndofs += press.assign_dofs(ndofs);

      sys.assemble();
      sys.solve(3, &xsln, &ysln, &psln);

     // visualization
      sprintf(title, "Velocity, time %g", current_time);
      vview.set_title(title);
      vview.show(&xsln, &ysln, EPS_LOW);
      sprintf(title, "Pressure, time %g", current_time);
      pview.set_title(title);
      pview.show(&psln);

      // solve fine problem
      RefSystem ref(&sys, 0); // just spacial refinement
      ref.assemble();
      ref.solve(3, &xref, &yref, &pref);

      // calculate errors
      DXDYFilter sln_vel(mag, &xsln, &ysln);
      DXDYFilter ref_vel(mag, &xref, &yref);

      double *crs_errors, *sln_errors;
      int    *crs_esort,  *sln_esort;

      double sln_err = 100 * calc_error(&sln_vel, &ref_vel, sln_esort, sln_errors);
      if (sln_err < space_tol || i == 1) done = true;
      info("Error %g%%", sln_err);


      if (done)
      {
        // solve super-coarse problem
        RefSystem crs(&sys, 0, -1);
        crs.assemble();
        crs.solve(3, &xcrs, &ycrs, &pcrs);

        DXDYFilter crs_vel(mag, &xcrs, &ycrs);
        double crs_err = 100 * calc_error(&crs_vel, &sln_vel, crs_esort, crs_errors);
      }

      // adapt the mesh (refine or coarse)
      adapt_mesh(done, &mesh, xcrs.get_mesh(), &xvel, crs_esort, crs_errors, sln_esort,  sln_errors);
      xvel.set_uniform_order(P_INIT_VEL);
      yvel.set_uniform_order(P_INIT_VEL);
      press.set_uniform_order(P_INIT_PRESSURE);

    }
    while(!done);


    xprev = xsln;
    yprev = ysln;
  }

  // wait for keyboard or mouse input
  View::wait("Waiting for all views to be closed.");
  return 0;
}
