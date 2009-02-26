#include "hermes2d.h"
#include "solver_umfpack.h"

// Reynolds numeber = L * U / vics, where
// L is real length of obstacle (in program always = 1)
// U is real input velocity (in program always = 1)
// visc is fluid viscosity
const double Re = 50000;

// initial polynomial orders
const int init_order_x = 2;
const int init_order_y = 2;
const int init_order_p = 1;

// time step size (dimensionless)
const double tau = 0.02;
// time when stop computation (dimensionless)
// real time = dimensionless time * (L / U)
const double final_time = 20.0;

// adaptivity parametes
const double delta = 0.1;   // tolerance for global error in percents
const double thr = 0.3;     // determine how many elements will be refine in one adaptive iteration
const double coef = 1.8;    // determine when to stop refinements within one adaptive iteration (error to be processed)

// stabilization
const bool stab = true;
const double delta_star = 1.0;
const double tau_star = 1.0;

// mesh file name
#define mesh_file "cylinder.mesh"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////   BOUNDARY CONDITIONS   ///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

const int BC_INLET = 1;
const int BC_OUTLET = 2;
const int BC_OUTER_TOP = 3;
const int BC_OUTER_BOTTOM = 4;
const int BC_WALL = 5;

int vel_bc_types(int marker)
  { return (marker != BC_OUTLET) ? BC_ESSENTIAL : BC_NONE; }

scalar xvel_bc_values(int marker, double x, double y)
  { return (marker != BC_WALL) ? 1 : 0; }

int press_bc_types(int marker)
  { return BC_NONE; }


////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////   WEAK FORMS   /////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

Solution xprev_crs, yprev_crs;
Solution xprev_sln, yprev_sln;
Solution xprev, yprev;

double *delta_K, *tau_K;

scalar bilinear_form_0_0_1_1_sym(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv) / Re +  
           int_u_v(fu, fv, ru, rv) / tau; }

scalar bilinear_form_0_0_1_1_unsym(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv); }
           
scalar bilinear_form_0_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return -int_u_dvdx(fu, fv, ru, rv); }

scalar bilinear_form_1_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return -int_u_dvdy(fu, fv, ru, rv); }

scalar linear_form_0(RealFunction* fv, RefMap* rv)
  { return int_u_v(&xprev, fv, rv, rv) / tau; }

scalar linear_form_1(RealFunction* fv, RefMap* rv)
  { return int_u_v(&yprev, fv, rv, rv) / tau; }


////////////////////   STABILIZATION   /////////////////////////////////////////////////////////////
#include "integrals_stab.cpp"

scalar bilinear_form_stab_0_0_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{ 
    double param = delta_K[ru->get_active_element()->id]; 
    return param * int_stab_0_0(&xprev, &yprev, fu, fv, ru, rv);
}

scalar bilinear_form_stab_0_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { 
    double param = tau_K[ru->get_active_element()->id];
    return param * int_dudx_dvdx(fu, fv, ru, rv); 
  }

scalar bilinear_form_stab_0_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { 
    double param = tau_K[ru->get_active_element()->id];
    return param * int_dudx_dvdy(fv, fu, rv, ru); 
  }

scalar bilinear_form_stab_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { 
    double param = tau_K[ru->get_active_element()->id];
    return param * int_dudy_dvdy(fu, fv, ru, rv); 
  }

scalar bilinear_form_stab_0_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { 
     double param = delta_K[ru->get_active_element()->id];
     return param * int_dudx_w_nabla_v(&xprev, &yprev, fu, fv, ru, rv); 
  }

scalar bilinear_form_stab_1_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { 
     double param = delta_K[ru->get_active_element()->id];
     return param * int_dudy_w_nabla_v(&xprev, &yprev, fu, fv, ru, rv); 
  }

scalar linear_form_stab_0(RealFunction* fv, RefMap* rv)
  { 
     RefMap* refmap = xprev.get_refmap();
     double param = delta_K[rv->get_active_element()->id];
     return param * int_w_nabla_u_v(&xprev, &yprev, fv, &xprev, rv, refmap) / tau;  
  }

scalar linear_form_stab_1(RealFunction* fv, RefMap* rv)
  { 
     RefMap* refmap = yprev.get_refmap();
     double param = delta_K[rv->get_active_element()->id];
     return param * int_w_nabla_u_v(&xprev, &yprev, fv, &yprev, rv, refmap) / tau;
  }

////////////  CALCULATION OF STAB PARAMETERS  //////////////////////////////////////////////////////
#include "stab_calculation.cpp"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////   ADAPTIVITY   ////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

static double* cmp_err;
static int compare(const void* p1, const void* p2)
{ 
  const int (*e1) = ((const int*) p1);
  const int (*e2) = ((const int*) p2);
  return cmp_err[(*e1)] < cmp_err[(*e2)] ? 1 : -1; 
}


double calc_error(MeshFunction* sln, MeshFunction* rsln, int*& esort, double*& errors)
{
  int i, j;

  int nact = sln->get_mesh()->get_num_active_elements();
  esort = new int[nact];

  double total_error = 0.0, total_norm = 0.0;
  int k = 0;

  Mesh* cmesh = sln->get_mesh();
  Mesh* fmesh = rsln->get_mesh();
    
  int max = cmesh->get_max_element_id();
  errors = new double[max];
  memset(errors, 0, sizeof(double) * max);
    
  Element* e;
  for_all_active_elements(e, cmesh)
  {
    sln->set_active_element(e);
    update_limit_table(e->get_mode());

    Element* fe = fmesh->get_element(e->id);
    if (fe->active) // both meshes are the same
    {
      rsln->set_active_element(fe);
      RefMap* crm = sln->get_refmap();
      RefMap* frm = rsln->get_refmap();

      double err  = int_l2_error(sln, rsln, crm, frm);
      total_norm += int_l2_norm (rsln, frm);

      errors[e->id] += err;
      total_error += err;
    }
    else // reference mesh IS refined
    {
      for (i = 0; i < 4; i++)
      {
        sln->push_transform(i);
        rsln->set_active_element(fe->sons[i]);
  
        RefMap* crm = sln->get_refmap();
        RefMap* frm = rsln->get_refmap();
  
        double err  = int_l2_error(sln, rsln, crm, frm);
        total_norm += int_l2_norm (rsln, frm);
  
        errors[e->id] += err;
        total_error += err;
  
        sln->pop_transform();
      }
    }       
    esort[k++] = e->id;
  }
  
  for_all_inactive_elements(e, cmesh)
    errors[e->id] = -1.0;

  for_all_active_elements(e, cmesh)
    errors[e->id] /= total_norm;

  cmp_err = errors;
  qsort(esort, nact, sizeof(int), compare);
    
  return sqrt(total_error / total_norm);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////////////////////////

void adapt_solution(bool& done, double to_be_processed, Mesh* mesh, Mesh* cmesh, Space* space, int* esort0, double* errors0, int* esort, double* errors)
{
  int i, j;
  double err0 = esort[0];
  double processed_error = 0.0;
  if (!done) // refinements
  {
    for (i = 0; i < mesh->get_num_active_elements(); i++)
    {
      int id = esort[i];
      double err = errors[id];
  
      if (err < thr * errors[esort[0]]) { break; }
      if ( (processed_error > coef * to_be_processed) && fabs((err - err0)/err0) > 1e-3 ) { info("Processed error condition."); break; }
      err0 = err;
      processed_error += err;
  
      Element* e;
      e = mesh->get_element(id);
      mesh->refine_element(id);
      for (j = 0; j < 4; j++)
        space->set_element_order(e->sons[j]->id, space->get_element_order(id));
    }
  }

  if (done)  // unrefinements
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

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////   MAIN PROGRAM   ////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  Mesh mesh, basemesh;
  basemesh.load(mesh_file);
  mesh.copy(&basemesh);
  
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  WeakForm wf(3);
  wf.add_biform(0, 0, bilinear_form_0_0_1_1_sym, SYM);
  wf.add_biform(0, 0, bilinear_form_0_0_1_1_unsym, UNSYM, ANY, 2, &xprev, &yprev);
  wf.add_biform(1, 1, bilinear_form_0_0_1_1_sym, SYM);
  wf.add_biform(1, 1, bilinear_form_0_0_1_1_unsym, UNSYM, ANY, 2, &xprev, &yprev);
  wf.add_biform(0, 2, bilinear_form_0_2, ANTISYM);
  wf.add_biform(1, 2, bilinear_form_1_2, ANTISYM);
  wf.add_liform(0, linear_form_0, ANY, 1, &xprev);
  wf.add_liform(1, linear_form_1, ANY, 1, &yprev);

  if (stab)
  {
    wf.add_biform(0, 0, bilinear_form_stab_0_0_1_1, UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(0, 0, bilinear_form_stab_0_0);
    wf.add_biform(0, 1, bilinear_form_stab_0_1, SYM);
    wf.add_biform(1, 1, bilinear_form_stab_0_0_1_1, UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(1, 1, bilinear_form_stab_1_1);
    wf.add_biform(0, 2, bilinear_form_stab_0_2, UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_biform(1, 2, bilinear_form_stab_1_2, UNSYM, ANY, 2, &xprev, &yprev);
    wf.add_liform(0, linear_form_stab_0, ANY, 2, &xprev, &yprev);
    wf.add_liform(1, linear_form_stab_1, ANY, 2, &xprev, &yprev);
  }
  
  H1Space xvel(&mesh, &shapeset);
  H1Space yvel(&mesh, &shapeset);
  H1Space press(&mesh, &shapeset);
  xvel.set_bc_types(vel_bc_types);
  xvel.set_bc_values(xvel_bc_values);
  yvel.set_bc_types(vel_bc_types);
  press.set_bc_types(press_bc_types);
  xvel.set_uniform_order(init_order_x);
  yvel.set_uniform_order(init_order_y);
  press.set_uniform_order(init_order_p);

  Mesh rmesh, cmesh;
  rmesh.copy(&mesh);  rmesh.refine_all_elements();
  cmesh.copy(&mesh);  cmesh.unrefine_all_elements();
  xprev_crs.set_const(&cmesh, 1.0);
  yprev_crs.set_zero(&cmesh);
  xprev_sln.set_const(&mesh, 1.0);
  yprev_sln.set_zero(&mesh);
  xprev.set_const(&rmesh, 1.0);
  yprev.set_zero(&rmesh);

  VectorView magview("Velocity", 0, 0, 800, 600);
  ScalarView pressview("Pressure", 900, 0, 800, 600);

  double total_time = 0.0;
  for (int it = 1; total_time < final_time; it++)
  {
    info("\n*** Time iteration %d, t = %g s ***\n", it, total_time += tau);

    Solution xcrs, ycrs, pcrs;
    Solution xsln, ysln, psln;
    Solution xref, yref, pref;

    bool done = false;
    int at = 0;
    do
    {
      info("\n*** Adaptive iteration %d ***\n", at++);

      rmesh.copy(&mesh);  rmesh.refine_all_elements();
      cmesh.copy(&mesh);  cmesh.unrefine_all_elements();

      // solve problem
      int ndofs = 0;
      ndofs += xvel.assign_dofs(ndofs);
      ndofs += yvel.assign_dofs(ndofs);
      ndofs += press.assign_dofs(ndofs);  
      if (stab) {
        int ne = mesh.get_max_element_id() + 1;
        delta_K = new double[ne];  tau_K = new double[ne];
        calculate_stabilization_parameters(delta_K, tau_K, &xprev_sln, &yprev_sln, &mesh);
      }
      UmfpackSolver umfpack;
      LinSystem sys(&wf, &umfpack);
      sys.set_spaces(3, &xvel, &yvel, &press);
      sys.set_pss(1, &pss);
      sys.assemble();
      sys.solve(3, &xsln, &ysln, &psln);
      if (stab) { delete [] delta_K;   delete [] tau_K; }
      
     // visualization
      magview.show(&xsln, &ysln, EPS_LOW);
      pressview.show(&psln);

      // solve fine problem
      if (stab) {
        int ne = rmesh.get_max_element_id() + 1;
        delta_K = new double[ne];  tau_K = new double[ne];
        calculate_stabilization_parameters(delta_K, tau_K, &xprev, &yprev, &rmesh);
      }
      RefSystem ref(&sys, 0);
      ref.assemble();
      ref.solve(3, &xref, &yref, &pref);
      if (stab) { delete [] delta_K;   delete [] tau_K; }

      // calculate errors
      DXDYFilter sln_vel(mag, &xsln, &ysln);
      DXDYFilter ref_vel(mag, &xref, &yref);

      double *crs_errors, *sln_errors;
      int    *crs_esort,  *sln_esort;

      double sln_err = 100 * calc_error(&sln_vel, &ref_vel, sln_esort, sln_errors);
      if (sln_err < delta) done = true;   
      info("Error %g%%", sln_err);    

      if (done) 
      {
        // solve super-coarse problem
        if (stab) {
          int ne = cmesh.get_max_element_id() + 1;
          delta_K = new double[ne];  tau_K = new double[ne];
          calculate_stabilization_parameters(delta_K, tau_K, &xprev_crs, &yprev_crs, &cmesh);
        }
        RefSystem crs(&sys, 0, -1);
        crs.assemble();
        crs.solve(3, &xcrs, &ycrs, &pcrs);
        if (stab) { delete [] delta_K;   delete [] tau_K; }

        DXDYFilter crs_vel(mag, &xcrs, &ycrs);
        double crs_err = 100 * calc_error(&crs_vel, &sln_vel, crs_esort, crs_errors);
      }

      // adapt the mesh
      double to_be_processed = sqr(sln_err/100) - sqr(delta/100);
      adapt_solution(done, to_be_processed, &mesh, xcrs.get_mesh(), &xvel, crs_esort, crs_errors, sln_esort,  sln_errors);
      xvel.set_uniform_order(init_order_x);
      yvel.set_uniform_order(init_order_y);
      press.set_uniform_order(init_order_p);
                  
    }
    while(!done);


    // update previous solutions
    xprev = xref;
    yprev = yref;

  }
  
  
  View::wait();
  return 0;
}
