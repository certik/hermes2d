
#include "hermes2d.h"
#include "solver_umfpack.h"

//// boundary conditions ///////////////////////////////////////////////////////////////////////////


const double Re = 1000000;
const double visc = 1.0/Re;
const double tau = 0.05;

// stabilization 
const bool stab = true;
const double delta_star = 1.0;
const double tau_star = 1.0;

// output:  marker = 3
// input:   marker = 2
// walls:   marker = 4
// airfoil: marker = 1

int xvel_bc_type(int marker)
  { return (marker != 3) ? BC_ESSENTIAL : BC_NONE; }

scalar xvel_bc_value(int marker, double x, double y)
  { return ((marker == 2) || (marker == 4)) ? 1 : 0; }

int yvel_bc_type(int marker)
  { return (marker != 3) ? BC_ESSENTIAL : BC_NONE; }

int press_bc_type(int marker)
  { return BC_NONE; }


////////////////////// Necessary integrals //////////////////////////////////////////

inline double int_stab_0_0(RealFunction* w1, RealFunction* w2,
                           RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = fu->get_fn_order() + fv->get_fn_order() + 2 * w1->get_fn_order() + ru->get_inv_ref_order(); 
  limit_order(o);

  w1->set_quad_order(o, FN_VAL);
  w2->set_quad_order(o, FN_VAL);
  fu->set_quad_order(o, FN_ALL);
  fv->set_quad_order(o);

  double *dvdx, *dvdy;
  fv->get_dx_dy_values(dvdx, dvdy);
  double* w1val = w1->get_fn_values();
  double* w2val = w2->get_fn_values();

  double *dudxx, *dudxy, *dudyy;
  dudxx = fu->get_dxx_values();
  dudxy = fu->get_dxy_values();    
  dudyy = fu->get_dyy_values();  
  double *dudx, *dudy;
  fu->get_dx_dy_values(dudx, dudy);
  double* uval = fu->get_fn_values();

  double result = 0.0;
  double3* pt = quad->get_points(o);
  int np = quad->get_num_points(o);
  double2x2 *mu, *mv;
  if (ru->is_jacobian_const()) 
  {
    for (int i = 0; i < np; i++) 
    {
      mu = ru->get_const_inv_ref_map(); 
      mv = rv->get_const_inv_ref_map(); 
      double a = (sqr((*mu)[0][0]) + sqr((*mu)[1][0]));
      double b = (sqr((*mu)[0][1]) + sqr((*mu)[1][1]));
      double c = 2.0 * ((*mu)[0][0]*(*mu)[0][1] + (*mu)[1][0]*(*mu)[1][1]);
      result += pt[i][2] * ( ( - ((dudxx[i]*a + dudxy[i]*c + dudyy[i]*b ) / Re)
                               + (w1val[i] * t_dudx + w2val[i] * t_dudy)
                               + (uval[i] / tau) 
                             )
                             * ((w1val[i] * t_dvdx + w2val[i] * t_dvdy)) ); 
    }
    result *= ru->get_const_jacobian();
  }
  else
  {
    mu = ru->get_inv_ref_map(o); 
    mv = rv->get_inv_ref_map(o); 
    double3x2 *mm;
    mm = ru->get_second_ref_map(o); 
    double* jac = ru->get_jacobian(o); 
    for (int i = 0; i < np; i++, mu++, mv++, mm++) 
    {
      double a = (sqr((*mu)[0][0]) + sqr((*mu)[1][0]));
      double b = (sqr((*mu)[0][1]) + sqr((*mu)[1][1]));
      double c = 2.0 * ((*mu)[0][0]*(*mu)[0][1] + (*mu)[1][0]*(*mu)[1][1]);
      double coefx = (*mm)[0][0] + (*mm)[2][0];
      double coefy = (*mm)[0][1] + (*mm)[2][1];
      result += pt[i][2] * (jac[i]) * ( ( - (( dudx[i]*coefx + dudy[i]*coefy + dudxx[i]*a + dudxy[i]*c + dudyy[i]*b ) / Re)
                                          + (w1val[i] * t_dudx + w2val[i] * t_dudy)
                                          + (uval[i] / tau) 
                                        )
                                        * ((w1val[i] * t_dvdx + w2val[i] * t_dvdy)) ); 
    }
  }
  return result;
}


inline double int_dudx_w_nabla_v(RealFunction* w1, RealFunction* w2,
                              RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = fu->get_fn_order() + fv->get_fn_order() + w1->get_fn_order() + ru->get_inv_ref_order(); 
  limit_order(o);

  w1->set_quad_order(o, FN_VAL);
  w2->set_quad_order(o, FN_VAL);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double *dudx, *dudy;
  fu->get_dx_dy_values(dudx, dudy);
  double *dvdx, *dvdy;
  fv->get_dx_dy_values(dvdx, dvdy);
  double* w1val = w1->get_fn_values();
  double* w2val = w2->get_fn_values();

  h1_integrate_dd_expression((t_dudx) * (w1val[i] * t_dvdx + w2val[i] * t_dvdy));
  return result;
}

inline double int_dudy_w_nabla_v(RealFunction* w1, RealFunction* w2,
                              RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  
  int o = fu->get_fn_order() + fv->get_fn_order() + w1->get_fn_order() + ru->get_inv_ref_order(); 
  limit_order(o);

  w1->set_quad_order(o, FN_VAL);
  w2->set_quad_order(o, FN_VAL);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  double *dudx, *dudy;
  fu->get_dx_dy_values(dudx, dudy);
  double *dvdx, *dvdy;
  fv->get_dx_dy_values(dvdx, dvdy);
  double* w1val = w1->get_fn_values();
  double* w2val = w2->get_fn_values();

  h1_integrate_dd_expression((t_dudy) * (w1val[i] * t_dvdx + w2val[i] * t_dvdy));
  return result;
}


inline double int_w_nabla_v_w(RealFunction* w1, RealFunction* w2,
                              RealFunction* fw, RealFunction* fv, RefMap* rw, RefMap* rv)
{
  Quad2D* quad = fv->get_quad_2d();

  int o = fw->get_fn_order() + fv->get_fn_order() + w1->get_fn_order() + rv->get_inv_ref_order(); 
  limit_order(o);

  w1->set_quad_order(o, FN_VAL);
  w2->set_quad_order(o, FN_VAL);
  fw->set_quad_order(o);
  fv->set_quad_order(o);

  double *dvdx, *dvdy;
  fv->get_dx_dy_values(dvdx, dvdy);
  double* wval = fw->get_fn_values();
  double* w1val = w1->get_fn_values();
  double* w2val = w2->get_fn_values();

  double result = 0.0;
  double3* pt = quad->get_points(o); 
  int np = quad->get_num_points(o);
  double2x2 *mv;
  mv = rv->get_inv_ref_map(o); 
  double* jac = rv->get_jacobian(o); 
  for (int i = 0; i < np; i++, mv++) 
    result += pt[i][2] * jac[i] * ((w1val[i] * t_dvdx + w2val[i] * t_dvdy) * wval[i]); 

  return result;
}

////////////////////////////////////////////////////////////////////////////////////

Solution xprev;
Solution yprev;

double *delta_K, *tau_K;

scalar bilinear_form_unsym_0_0_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv) / Re +
           int_u_v(fu, fv, ru, rv) / tau +
           int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv); }

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

scalar bilinear_form_stab_1_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { 
    double param = tau_K[ru->get_active_element()->id];
    return param * int_dudx_dvdy(fu, fv, ru, rv); 
  }

scalar bilinear_form_stab_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { 
    double param = tau_K[ru->get_active_element()->id];
    return param * int_dudy_dvdy(fu, fv, ru, rv); 
  }

scalar bilinear_form_unsym_0_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { 
     if (stab) {
       double param = delta_K[ru->get_active_element()->id];
       return -int_u_dvdx(fu, fv, ru, rv) + param * int_dudx_w_nabla_v(&xprev, &yprev, fu, fv, ru, rv); 
     }
     else return -int_u_dvdx(fu, fv, ru, rv);
  }

scalar bilinear_form_unsym_1_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { 
     if (stab) {
       double param = delta_K[ru->get_active_element()->id];
       return -int_u_dvdy(fu, fv, ru, rv) + param * int_dudy_w_nabla_v(&xprev, &yprev, fu, fv, ru, rv); 
     }
     else return -int_u_dvdy(fu, fv, ru, rv);
  }

scalar bilinear_form_unsym_2_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_dudx_v(fu, fv, ru, rv); }

scalar bilinear_form_unsym_2_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_dudy_v(fu, fv, ru, rv); }

scalar linear_form_0(RealFunction* fv, RefMap* rv)
  { 
     RefMap* refmap = xprev.get_refmap();
     if (stab) {
       double param = delta_K[rv->get_active_element()->id];
       return int_u_v(&xprev, fv, refmap, rv) / tau + param * int_w_nabla_v_w(&xprev, &yprev, &xprev, fv, refmap, rv) / tau;  
     }
     else return int_u_v(&xprev, fv, refmap, rv) / tau; 
  }

scalar linear_form_1(RealFunction* fv, RefMap* rv)
  { 
     RefMap* refmap = yprev.get_refmap();
     if (stab) {
       double param = delta_K[rv->get_active_element()->id];
       return int_u_v(&yprev, fv, refmap, rv) / tau + param * int_w_nabla_v_w(&xprev, &yprev, &yprev, fv, refmap, rv) / tau;
     }
     else return int_u_v(&yprev, fv, refmap, rv) / tau;
  }



////////////////////////////////////////////////////////////////////////////////////////////////////
//// Calculation of stabilization parameters //////////////////////////////////////////////////////

void calculate_elements_length(double* ele_len, double* u_infty, Solution* u1, Solution* u2, Mesh* mesh)
{
  Quad2D* quad = &g_quad_2d_std;
  u1->set_quad_2d(quad);
  u2->set_quad_2d(quad);
  
  Solution tmp;
  tmp.set_zero(mesh);

  Mesh* meshes[3] = { mesh, u1->get_mesh(), u2->get_mesh() };
  Transformable* tr[3] = { &tmp, u1, u2 };
  Traverse trav;
  trav.begin(3, meshes, tr);
  
  int ne = mesh->get_max_element_id() + 1;
  double* a = new double[ne];
  double* b = new double[ne];
  int* n = new int[ne];
  memset(a, 0, sizeof(double) * ne);
  memset(b, 0, sizeof(double) * ne);
  memset(n, 0, sizeof(int) * ne);

  Element** ee; 
  while ((ee = trav.get_next_state(NULL, NULL)) != NULL)
  {
    int o = u1->get_fn_order();
    u1->set_quad_order(o);
    u2->set_quad_order(o);
    scalar *uval, *vval;
    uval = u1->get_fn_values();
    vval = u2->get_fn_values();
    int np = quad->get_num_points(o);
    for (int i = 0; i < np; i++) {
      a[ee[0]->id] += uval[i];
      b[ee[0]->id] += vval[i];
    }
    n[ee[0]->id] += np;
  }
  trav.finish();

  Element* e;
  for_all_active_elements(e, mesh)
  {
    // averaged values of velocities on elements of mesh
    double c = a[e->id] / n[e->id]; 
    double d = b[e->id] / n[e->id]; 
    u_infty[e->id] = (fabs(c) > fabs(d)) ? fabs(c) : fabs(d);

    double den = (sqrt(sqr(c) + sqr(d)));
    // find the max length
    double min = 10000.0, max = -10000.0;
    for (int i = 0; i < e->nvert; i++)
    {
      double x = e->vn[i]->x;
      double y = e->vn[i]->y;
      double l = (c*x + d*y) / den;
      if (l > max) max = l;
      if (l < min) min = l;
    }

    ele_len[e->id] = max - min;
  }

  delete [] a;
  delete [] b;
  delete [] n;

}


void calculate_stabilization_parameters(double* delta_K, double* tau_K, double* h, double* u_infty, Mesh* mesh)
{
  Element* e;
  for_all_active_elements(e, mesh)
  {
    delta_K[e->id] = delta_star * sqr(h[e->id]);
    tau_K[e->id] = tau_star * 1.0;   
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// Refinement criterions //////////////////////////////////////////////////////////////

int criterion1(Element* e)
{
  if (e->is_triangle()) return -1;

  double mid1x, mid1y, mid2x, mid2y;
  mid1x = 0.5 * (e->vn[1]->x + e->vn[2]->x);
  mid1y = 0.5 * (e->vn[1]->y + e->vn[2]->y);
  mid2x = 0.5 * (e->vn[3]->x + e->vn[0]->x);
  mid2y = 0.5 * (e->vn[3]->y + e->vn[0]->y);
  double h = sqrt(sqr(mid1x - mid2x) + sqr(mid1y - mid2y));
  mid1x = 0.5 * (e->vn[0]->x + e->vn[1]->x);
  mid1y = 0.5 * (e->vn[0]->y + e->vn[1]->y);
  mid2x = 0.5 * (e->vn[2]->x + e->vn[3]->x);
  mid2y = 0.5 * (e->vn[2]->y + e->vn[3]->y);
  double v = sqrt(sqr(mid1x - mid2x) + sqr(mid1y - mid2y));

  if (h > 3 * v) return 2;
  if (v > 3 * h) return 1;

  return -1;
}


int criterion2(Element* e)
{
  for (int i = 0; i < e->nvert; i++)
    if ((e->vn[i]->x < 1.5) && (e->vn[i]->x > -0.2) && (e->vn[i]->y < 0.5) && (e->vn[i]->y > -0.5)) 
    {  
      if (e->is_triangle()) return 0;
    
      double mid1x, mid1y, mid2x, mid2y;
      mid1x = 0.5 * (e->vn[1]->x + e->vn[2]->x);
      mid1y = 0.5 * (e->vn[1]->y + e->vn[2]->y);
      mid2x = 0.5 * (e->vn[3]->x + e->vn[0]->x);
      mid2y = 0.5 * (e->vn[3]->y + e->vn[0]->y);
      double h = sqrt(sqr(mid1x - mid2x) + sqr(mid1y - mid2y));
      mid1x = 0.5 * (e->vn[0]->x + e->vn[1]->x);
      mid1y = 0.5 * (e->vn[0]->y + e->vn[1]->y);
      mid2x = 0.5 * (e->vn[2]->x + e->vn[3]->x);
      mid2y = 0.5 * (e->vn[2]->y + e->vn[3]->y);
      double v = sqrt(sqr(mid1x - mid2x) + sqr(mid1y - mid2y));
    
      if (h > 3 * v) return 2;
      if (v > 3 * h) return 1;
      return 0;
    }

  return -1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);
  
  H1ShapesetBeuchler shapeset;
  PrecalcShapeset pss(&shapeset);
  
  Mesh mesh;
  mesh.load("airfoil.mesh");

  mesh.refine_by_criterion(criterion1, 2);  
  mesh.refine_by_criterion(criterion2, 3);
  mesh.refine_towards_boundary(1, 1);

  H1Space xvel(&mesh, &shapeset);
  xvel.set_bc_types(xvel_bc_type);
  xvel.set_bc_values(xvel_bc_value);
  xvel.set_uniform_order(2);
  
  H1Space yvel(&mesh, &shapeset);
  yvel.set_bc_types(yvel_bc_type);
  yvel.set_uniform_order(2);
  
  H1Space press(&mesh, &shapeset);
  press.set_bc_types(press_bc_type);
  press.set_uniform_order(1);

  xprev.set_const(&mesh, 1.0);
  yprev.set_const(&mesh, 0.0);

 
  WeakForm wf(3);
  wf.add_biform(0, 0, bilinear_form_unsym_0_0_1_1, UNSYM, 0, 2, &xprev, &yprev);
  if (stab) {
    wf.add_biform(0, 0, bilinear_form_stab_0_0_1_1, UNSYM, 0, 2, &xprev, &yprev);
    wf.add_biform(0, 0, bilinear_form_stab_0_0);
    wf.add_biform(0, 1, bilinear_form_stab_0_1);
  }
  wf.add_biform(0, 2, bilinear_form_unsym_0_2, UNSYM, 0, 2, &xprev, &yprev);

  wf.add_biform(1, 1, bilinear_form_unsym_0_0_1_1, UNSYM, 0, 2, &xprev, &yprev);
  if (stab) {
    wf.add_biform(1, 1, bilinear_form_stab_0_0_1_1, UNSYM, 0, 2, &xprev, &yprev);
    wf.add_biform(1, 1, bilinear_form_stab_1_1);
    wf.add_biform(1, 0, bilinear_form_stab_1_0);
  }
  wf.add_biform(1, 2, bilinear_form_unsym_1_2, UNSYM, 0, 2, &xprev, &yprev);

  wf.add_biform(2, 0, bilinear_form_unsym_2_0);
  wf.add_biform(2, 1, bilinear_form_unsym_2_1);
  wf.add_liform(0, linear_form_0, 0, 2, &xprev, &yprev);
  wf.add_liform(1, linear_form_1, 0, 2, &xprev, &yprev);


  // set up solver and linear system
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);

  Solution xsln, ysln, psln;

  VectorView vview("Velocity", 0, 0, 600, 500);
  ScalarView pview("Pressure", 600, 0, 600, 500);

  for (int i = 0; i < 1000; i++)
  {
    printf("\n-------------- Iteration %d, t = %g s -----------------------\n", i, i*tau);

    // calculate stabilization parameters
    int ne = mesh.get_max_element_id() + 1;
    double* ele_len = new double[ne];
    double* u_infty = new double[ne];
    delta_K = new double[ne];
    tau_K = new double[ne];
    if (stab) {
      calculate_elements_length(ele_len, u_infty, &xprev, &yprev, &mesh); 
      calculate_stabilization_parameters(delta_K, tau_K, ele_len, u_infty, &mesh);
    }
    delete [] ele_len;
    delete [] u_infty;

    int ndofs = 0;
    ndofs += xvel.assign_dofs(ndofs);
    ndofs += yvel.assign_dofs(ndofs);
    ndofs += press.assign_dofs(ndofs);

    sys.set_spaces(3, &xvel, &yvel, &press);
    sys.set_pss(1, &pss);
    sys.assemble();
    sys.solve(3, &xsln, &ysln, &psln);

    delete [] delta_K;
    delete [] tau_K;

    pview.set_min_max_range(-0.5,0.5);
    pview.show_contours(0.02);
    pview.show(&psln);
    vview.set_min_max_range(0.0,1.2);
    vview.show(&xsln, &ysln);


    xprev = xsln;
    yprev = ysln;

    printf("\n");
  }
  
  hermes2d_finalize();
}

