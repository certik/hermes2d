//#define DEBUG_ORDER
#include "hermes2d.h"


const double tau = 0.01;

const double gam = 1.4;

const double rho_left = 1.0;
const double rho_top = 1.7;

const double v1_left = 2.9;
const double v2_left = 0.0;
const double v1_top = 2.619334;
const double v2_top = -0.5063;

const double p_left = 0.714286;
const double p_top = 1.52819;


//// boundary conditions ///////////////////////////////////////////////////////////////////////////

double zhlazeni1(double a, double b, double x)
{
  return 0.5 * (a+b) + 0.5 * (a-b) * (2.0 / M_PI) * atan(1000 * (1-x));  
}

double zhlazeni2(double a, double b, double x)
{
  return 0.5 * (a+b) + 0.5 * (b-a) * (2.0 / M_PI) * atan(1000 * (x));  
}

// BC for density
int dens_bc_types(int marker)
{
  if (marker == 1 || marker == 2) return BC_ESSENTIAL;
  else if (marker == 3) return BC_NONE;
  else if (marker == 4) return BC_NATURAL;
}

scalar dens_bc_values(int marker, double x, double y) 
{
  //if (marker == 1) return rho_left;
  if (marker == 1) return zhlazeni1(rho_left, rho_top, y);
  //else if (marker == 2) return rho_top;
  else if (marker == 2) return zhlazeni2(rho_left, rho_top, x);
}

// BC for momentum #1
int mom1_bc_types(int marker)
{
  if (marker == 1 || marker == 2) return BC_ESSENTIAL;
  else if (marker == 3) return BC_NONE;
  else if (marker == 4) return BC_NATURAL;

}

scalar mom1_bc_values(int marker, double x, double y)
{
/*  if (marker == 1) return rho_left * v1_left;
  else if (marker == 2) return rho_top * v1_top;*/
  if (marker == 1) return zhlazeni1(rho_left * v1_left, rho_top * v1_top, y);
  else if (marker == 2) return zhlazeni2(rho_left * v1_left, rho_top * v1_top, x);

}

// BC for momentum #2
int mom2_bc_types(int marker)
{
  if (marker == 1 || marker == 2 || marker == 3) return BC_ESSENTIAL;
  else if (marker == 4) return BC_NATURAL;
}

scalar mom2_bc_values(int marker, double x, double y)
{
/*  if (marker == 1) return v2_left * rho_left;
  else if (marker == 2) return v2_top * rho_top;*/
  if (marker == 1) return zhlazeni1(rho_left * v2_left, rho_top * v2_top, y);
  else if (marker == 2) return zhlazeni2(rho_left * v2_left, rho_top * v2_top, x);
  else if (marker == 3) return 0.0;
}

// BC for energy
int energ_bc_types(int marker)
{
  if (marker == 1 || marker == 2) return BC_ESSENTIAL;
  else if (marker == 3) return BC_NONE;
  else if (marker == 4) return BC_NATURAL;
}

scalar energ_bc_values(int marker, double x, double y)
{
/*  if (marker == 1) return (p_left / (gam - 1)) + (rho_left *(v1_left*v1_left + v2_left*v2_left)) / 2;
  else if (marker == 2) return (p_top / (gam - 1)) + (rho_top *(v1_top*v1_top + v2_top*v2_top)) / 2;*/
  if (marker == 1) return zhlazeni1((p_left / (gam - 1)) + (rho_left *(v1_left*v1_left + v2_left*v2_left)) / 2, 
                                    (p_top / (gam - 1)) + (rho_top *(v1_top*v1_top + v2_top*v2_top)) / 2, y);
  else if (marker == 2) return zhlazeni2((p_left / (gam - 1)) + (rho_left *(v1_left*v1_left + v2_left*v2_left)) / 2, 
                                         (p_top / (gam - 1)) + (rho_top *(v1_top*v1_top + v2_top*v2_top)) / 2, x);
}




//// integrals /////////////////////////////////////////////////////////////////////////////////////

inline double int_w_dudx_v(RealFunction* w, RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = /*w->get_fn_order()*/ 2 + fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  w->set_quad_order(o, FN_VAL);
  fu->set_quad_order(o);
  fv->set_quad_order(o, FN_VAL);

  double *wval = w->get_fn_values();  
  double *dudx = fu->get_dx_values();
  double *dudy = fu->get_dy_values();
  double *vval = fv->get_fn_values();

  h1_integrate_dd_expression(wval[i] * t_dudx * vval[i]);
  return result;
}

inline double int_w_dudy_v(RealFunction* w, RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = /*w->get_fn_order()*/ 2 + fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  w->set_quad_order(o, FN_VAL);
  fu->set_quad_order(o);
  fv->set_quad_order(o, FN_VAL);

  double *wval = w->get_fn_values();  
  double *dudx = fu->get_dx_values();
  double *dudy = fu->get_dy_values();
  double *vval = fv->get_fn_values();

  h1_integrate_dd_expression(wval[i] * t_dudy * vval[i]);
  return result;
}

  
//// bilinear and linear forms /////////////////////////////////////////////////////////////////////

Solution u1prev, u2prev, u3prev, u4prev;

SimpleFilter *f10a, *f10b, *f11a, *f11b, *f12a, *f12b;
SimpleFilter *f20a, *f20b, *f21b, *f22b;
SimpleFilter *f30a, *f30b, *f31a, *f31b, *f32b, *f33a, *f33b;

// ROW 0
scalar bilinear_form_0_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv) 
  { return int_u_v(fu, fv, ru, rv) / tau; }

scalar bilinear_form_0_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv) 
  { return int_dudx_v(fu, fv, ru, rv); }

scalar bilinear_form_0_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_dudy_v(fu, fv, ru, rv); }

// ROW 1
scalar bilinear_form_1_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_w_dudx_v(f10a, fu, fv, ru, rv) + int_w_dudy_v(f10b, fu, fv, ru, rv); }

scalar bilinear_form_1_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_u_v(fu, fv, ru, rv) / tau + int_w_dudx_v(f11a, fu, fv, ru, rv) + int_w_dudy_v(f11b, fu, fv, ru, rv); }

scalar bilinear_form_1_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_w_dudx_v(f12a, fu, fv, ru, rv) + int_w_dudy_v(f12b, fu, fv, ru, rv); }

scalar bilinear_form_1_3(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return (gam - 1) * int_dudx_v(fu, fv, ru, rv); }
  
// ROW 2
scalar bilinear_form_2_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv) 
  { return int_w_dudx_v(f20a, fu, fv, ru, rv) + int_w_dudy_v(f20b, fu, fv, ru, rv); }

scalar bilinear_form_2_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv) 
  { return int_w_dudx_v(f11b, fu, fv, ru, rv) + int_w_dudy_v(f21b, fu, fv, ru, rv); }

scalar bilinear_form_2_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_u_v(fu, fv, ru, rv) / tau + int_w_dudx_v(f12b, fu, fv, ru, rv) + int_w_dudy_v(f22b, fu, fv, ru, rv); }

scalar bilinear_form_2_3(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return (gam - 1) * int_dudy_v(fu, fv, ru, rv); }
  
// ROW 3
scalar bilinear_form_3_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_w_dudx_v(f30a, fu, fv, ru, rv) + int_w_dudy_v(f30b, fu, fv, ru, rv); }

scalar bilinear_form_3_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_w_dudx_v(f31a, fu, fv, ru, rv) + int_w_dudy_v(f31b, fu, fv, ru, rv); }

scalar bilinear_form_3_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_w_dudx_v(f31b, fu, fv, ru, rv) + int_w_dudy_v(f32b, fu, fv, ru, rv); }

scalar bilinear_form_3_3(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv) 
  { return int_u_v(fu, fv, ru, rv) / tau + int_w_dudx_v(f33a, fu, fv, ru, rv) + int_w_dudy_v(f33b, fu, fv, ru, rv); }

    
//// linear forms //////////////////////////////////////////////////////////////////////////////////

scalar linear_form_0(RealFunction* fv, RefMap* rv) 
  { return int_u_v(&u1prev, fv, rv, rv) / tau; }

scalar linear_form_1(RealFunction* fv, RefMap* rv)
  { return int_u_v(&u2prev, fv, rv, rv) / tau; }

scalar linear_form_2(RealFunction* fv, RefMap* rv)
  { return int_u_v(&u3prev, fv, rv, rv) / tau; }

scalar linear_form_3(RealFunction* fv, RefMap* rv)
  { return int_u_v(&u4prev, fv, rv, rv) / tau; }


//// filters ///////////////////////////////////////////////////////////////////////////////////////

#define def_flt_2(name, var1, var2, exp) \
  void name(int n, scalar* var1, scalar* var2, scalar* rslt) \
    { for (int i = 0; i < n; i++) rslt[i] = exp; }
    
#define def_flt_3(name, var1, var2, var3, exp) \
  void name(int n, scalar* var1, scalar* var2, scalar* var3, scalar* rslt) \
    { for (int i = 0; i < n; i++) rslt[i] = exp; }
    
#define def_flt_4(name, var1, var2, var3, var4, exp) \
  void name(int n, scalar* var1, scalar* var2, scalar* var3, scalar* var4, scalar* rslt) \
    { for (int i = 0; i < n; i++) rslt[i] = exp; }

    
def_flt_3(fn10a, u1, u2, u3,  (gam - 1) * (0.5 *(sqr(u2[i]) + sqr(u3[i])) / sqr(u1[i])) - sqr(u2[i])/sqr(u1[i]));
def_flt_3(fn10b, u1, u2, u3, -(u2[i] * u3[i]) / sqr(u1[i]));
def_flt_2(fn11a, u1, u2,      (3 - gam) * u2[i] / u1[i]);
def_flt_2(fn11b, u1, u3,      u3[i] / u1[i]);
def_flt_2(fn12a, u1, u3,      (1 - gam) * u3[i] / u1[i]);
def_flt_2(fn12b, u1, u2,      u2[i] / u1[i]);

def_flt_3(fn20a, u1, u2, u3, -(u2[i] * u3[i]) / sqr(u1[i]));
def_flt_3(fn20b, u1, u2, u3,  (gam - 1) * (0.5 *(sqr(u2[i]) + sqr(u3[i])) / sqr(u1[i])) - sqr(u3[i])/sqr(u1[i]));
def_flt_2(fn21b, u1, u2,      (1 - gam) * u2[i] / u1[i]);
def_flt_2(fn22b, u1, u3,      (3 - gam) * u3[i] / u1[i]);    

def_flt_4(fn30a, u1, u2, u3, u4,  (u2[i] / u1[i]) * ((gam - 1) * ((sqr(u2[i]) + sqr(u3[i])) / sqr(u1[i])) - gam * u4[i]/u1[i]));
def_flt_4(fn30b, u1, u2, u3, u4,  (u3[i] / u1[i]) * ((gam - 1) * ((sqr(u2[i]) + sqr(u3[i])) / sqr(u1[i])) - gam * u4[i]/u1[i]));
def_flt_4(fn31a, u1, u2, u3, u4,  gam * u4[i]/u1[i] - (gam -1) * (sqr(u2[i])/sqr(u1[i])) - 0.5 * (gam - 1) * ((sqr(u2[i]) + sqr(u3[i])) / sqr(u1[i])));
def_flt_3(fn31b, u1, u2, u3,      (1 - gam) * (u2[i] * u3[i] / sqr(u1[i])));
def_flt_4(fn32b, u1, u2, u3, u4,  gam * u4[i]/u1[i] - (gam -1) * (sqr(u3[i])/sqr(u1[i])) - 0.5 * (gam - 1) * ((sqr(u2[i]) + sqr(u3[i])) / sqr(u1[i])));
def_flt_2(fn33a, u1, u2,          gam * u2[i]/u1[i]);
def_flt_2(fn33b, u1, u3,          gam * u3[i]/u1[i]);

////////////////////////////////////////////////////////////////////////////////////////////////////

int criterion(Element* e)
{

  int region[4];
  for (int i = 0; i < e->nvert; i++)
  {
    if (e->vn[i]->y < -0.554 * e->vn[i]->x + 0.9  ) region[i] = 1;
    else if (e->vn[i]->y <  0.435 * e->vn[i]->x - 0.886) region[i] = 3;
    else if ((e->vn[i]->y >  0.435 * e->vn[i]->x - 0.686) && (e->vn[i]->y > -0.554 * e->vn[i]->x + 1.1)) region[i] = 2;
    else region[i] = 0;    
  }
  for(int i = 0; i < e->nvert; i++)
    if (region[i] == 0) return 0;

  return -1;

}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("sq.mesh");
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  mesh.refine_all_elements();
  //mesh.refine_all_elements();
  //mesh.refine_all_elements();
  //mesh.refine_all_elements();
  //mesh.refine_by_criterion(criterion, 3);
  //mesh.regularize(0);
  

  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);

  // spaces
  H1Space space1(&mesh, &shapeset);
  H1Space space2(&mesh, &shapeset);
  H1Space space3(&mesh, &shapeset);
  H1Space space4(&mesh, &shapeset);

  // initialize boundary conditions
  space1.set_bc_types(dens_bc_types);
  space1.set_bc_values(dens_bc_values);
  space2.set_bc_types(mom1_bc_types);
  space2.set_bc_values(mom1_bc_values);
  space3.set_bc_types(mom2_bc_types);
  space3.set_bc_values(mom2_bc_values);
  space4.set_bc_types(energ_bc_types);
  space4.set_bc_values(energ_bc_values);

  // set polynomial degrees
  space1.set_uniform_order(1);
  space2.set_uniform_order(1);
  space3.set_uniform_order(1);
  space4.set_uniform_order(1);

  // assign degrees of freedom
  int ndofs = 0;
  ndofs += space1.assign_dofs(ndofs);
  ndofs += space2.assign_dofs(ndofs);
  ndofs += space3.assign_dofs(ndofs);
  ndofs += space4.assign_dofs(ndofs);

  // set initial values
  u1prev.set_const(&mesh, 1.0);
  u2prev.set_const(&mesh, 2.9);
  u3prev.set_zero(&mesh);
  u4prev.set_const(&mesh, 5.99);

  // problem initalization
  DiscreteProblem dp;
  dp.set_num_equations(4);
  dp.set_spaces(4, &space1, &space2, &space3, &space4);
  dp.set_pss(1, &pss);

  // set up weak formulation
  dp.set_bilinear_form(0, 0, bilinear_form_0_0);
  dp.set_bilinear_form(0, 1, bilinear_form_0_1);
  dp.set_bilinear_form(0, 2, bilinear_form_0_2);
  dp.set_bilinear_form(1, 0, bilinear_form_1_0);
  dp.set_bilinear_form(1, 1, bilinear_form_1_1);
  dp.set_bilinear_form(1, 2, bilinear_form_1_2);
  dp.set_bilinear_form(1, 3, bilinear_form_1_3);
  dp.set_bilinear_form(2, 0, bilinear_form_2_0);
  dp.set_bilinear_form(2, 1, bilinear_form_2_1);
  dp.set_bilinear_form(2, 2, bilinear_form_2_2);
  dp.set_bilinear_form(2, 3, bilinear_form_2_3);
  dp.set_bilinear_form(3, 0, bilinear_form_3_0);
  dp.set_bilinear_form(3, 1, bilinear_form_3_1);
  dp.set_bilinear_form(3, 2, bilinear_form_3_2);
  dp.set_bilinear_form(3, 3, bilinear_form_3_3);

  dp.set_linear_form(0, linear_form_0);
  dp.set_linear_form(1, linear_form_1);
  dp.set_linear_form(2, linear_form_2);
  dp.set_linear_form(3, linear_form_3);

  // visualization
  ScalarView dview("Density", 25, 0, 1550, 380);
  VectorView vview("Velocity", 25, 500, 1550, 390);

  // precalculate matrix sparse structure, allocate
  dp.create_matrix();

  // main loop
  for (int i = 0; i < 1000; i++)
  {
    printf("\n*** Iteration %d ***\n", i);

    // set up filters
    SimpleFilter filter10a(fn10a, &u1prev, &u2prev, &u3prev);
    SimpleFilter filter10b(fn10b, &u1prev, &u2prev, &u3prev);
    SimpleFilter filter11a(fn11a, &u1prev, &u2prev);
    SimpleFilter filter11b(fn11b, &u1prev, &u3prev);
    SimpleFilter filter12a(fn12a, &u1prev, &u3prev);
    SimpleFilter filter12b(fn12b, &u1prev, &u2prev);

    SimpleFilter filter20a(fn20a, &u1prev, &u2prev, &u3prev);
    SimpleFilter filter20b(fn20b, &u1prev, &u2prev, &u3prev);
    SimpleFilter filter21b(fn21b, &u1prev, &u2prev);
    SimpleFilter filter22b(fn22b, &u1prev, &u3prev);

    SimpleFilter filter30a(fn30a, &u1prev, &u2prev, &u3prev, &u4prev);
    SimpleFilter filter30b(fn30b, &u1prev, &u2prev, &u3prev, &u4prev);
    SimpleFilter filter31a(fn31a, &u1prev, &u2prev, &u3prev, &u4prev);
    SimpleFilter filter31b(fn31b, &u1prev, &u2prev, &u3prev);
    SimpleFilter filter32b(fn32b, &u1prev, &u2prev, &u3prev, &u4prev);
    SimpleFilter filter33a(fn33a, &u1prev, &u2prev);
    SimpleFilter filter33b(fn33b, &u1prev, &u3prev);

    // visualization
    dview.set_min_max_range(0.9,2.8);
    dview.show(&u1prev); 
    //vview.show(&filter12b, &filter11b);
    
    // assemble and solve
    dp.set_external_fns(21, &u1prev, &u2prev, &u3prev, &u4prev,
                            f10a = &filter10a,  f10b = &filter10b,
                            f11a = &filter11a,  f11b = &filter11b,
                            f12a = &filter12a,  f12b = &filter12b,
                            f20a = &filter20a,  f20b = &filter20b,
                            f21b = &filter21b,  f22b = &filter22b,
                            f30a = &filter30a,  f30b = &filter30b,
                            f31a = &filter31a,  f31b = &filter31b,
                                                f32b = &filter32b,
                            f33a = &filter33a,  f33b = &filter33b);

    dp.assemble_matrix_and_rhs();
    dp.solve_system(4, &u1prev, &u2prev, &u3prev, &u4prev);

  }

  // done
  //View::wait();
}

