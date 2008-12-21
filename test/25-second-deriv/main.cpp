// transformation of second derivatives in weak formulation
// laplace equation on one element with nonconstant jacobian solved twice 
// - standard weak form    ( grad u , grad v ) = v
// - using 2nd derivatives ( laplace u , v )   = v

#include "hermes2d.h"
#include "solver_umfpack.h"


inline double int_laplace_u_v(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv, bool adapt)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order(); 
  limit_order(o);

  fu->set_quad_order(o, FN_ALL);
  fv->set_quad_order(o, FN_VAL);

  double *vval = fv->get_fn_values();
  double *dudxx, *dudxy, *dudyy;
  dudxx = fu->get_dxx_values();
  dudxy = fu->get_dxy_values();    
  dudyy = fu->get_dyy_values();  
  double *dudx, *dudy;
  fu->get_dx_dy_values(dudx, dudy);

  double result = 0.0;
  double3* pt = quad->get_points(o);
  int np = quad->get_num_points(o);
  double2x2 *mu;
  mu = ru->get_inv_ref_map(o); 
  if (ru->is_jacobian_const()) 
  {
    for (int i = 0; i < np; i++, mu++) 
    {
      double a = (sqr((*mu)[0][0]) + sqr((*mu)[1][0]));
      double b = (sqr((*mu)[0][1]) + sqr((*mu)[1][1]));
      double c = 2.0 * ((*mu)[0][0]*(*mu)[0][1] + (*mu)[1][0]*(*mu)[1][1]);
      result += pt[i][2] * (dudxx[i]*a + dudxy[i]*c + dudyy[i]*b ) * vval[i]; 
    }
    result *= ru->get_const_jacobian();
  }
  else
  {
    double3x2 *mm;
    mm = ru->get_second_ref_map(o); 
    double* jac = ru->get_jacobian(o); 
    for (int i = 0; i < np; i++, mu++, mm++) 
    {
      double a = (sqr((*mu)[0][0]) + sqr((*mu)[1][0]));
      double b = (sqr((*mu)[0][1]) + sqr((*mu)[1][1]));
      double c = 2.0 * ((*mu)[0][0]*(*mu)[0][1] + (*mu)[1][0]*(*mu)[1][1]);
      double coefx = (*mm)[0][0] + (*mm)[2][0];
      double coefy = (*mm)[0][1] + (*mm)[2][1];
      result += pt[i][2] * (jac[i]) * (dudx[i]*coefx + dudy[i]*coefy + dudxx[i]*a + dudxy[i]*c + dudyy[i]*b) * vval[i] ; 
    }
  }


  if (adapt && (fabs(result) > 1e-5))
  {
    scalar sum = 0.0;
    for (int i = 0; i < 4; i++)
    {
      fu->push_transform(i);
      ru->push_transform(i);
      fv->push_transform(i);
      sum += int_laplace_u_v(fu, fv, ru, rv, false);
      fu->pop_transform();
      ru->pop_transform();
      fv->pop_transform();
    } 
    
    if (fabs(sum - result) / fabs(sum) < 1e-3) return sum;
    
    if (fv->get_ctm()->m[0] < 1.0 / 256.0)
    {
      warn("Adaptive quadrature: could not reach required accuracy.");
      return sum;
    }
    
    sum = 0.0;
    for (int i = 0; i < 4; i++)
    {
      fu->push_transform(i);
      ru->push_transform(i);
      fv->push_transform(i);
      sum += int_laplace_u_v(fu, fv, ru, rv, true);
      fu->pop_transform();
      ru->pop_transform();
      fv->pop_transform();
    }
    
    return sum;
  }
  
  return result;
}


scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_grad_u_grad_v(fu, fv, ru, rv);
}


scalar bilinear_form2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return -int_laplace_u_v(fu, fv, ru, rv, true);
}

scalar linear_form(RealFunction* fv, RefMap* rv)
{
  return int_v(fv, rv);
}


int main(int argc, char* argv[])
{
  Mesh mesh;
  H1ShapesetBeuchler shapeset;
  PrecalcShapeset pss(&shapeset);
  H1Space space(&mesh, &shapeset);
  Solution sln, sln2;
  UmfpackSolver umfpack;

  WeakForm wf(1);
  wf.add_biform(0, 0, bilinear_form);
  wf.add_liform(0, linear_form);

  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);

  mesh.load("square1.mesh");

  space.set_uniform_order(5);
  space.assign_dofs();

  sys.assemble();
  sys.solve(1, &sln);

  ScalarView view1("Solution 1 - grad grad");
  view1.show(&sln);

////////////////////////////////////////////////

  mesh.load("square1.mesh");
  
  Mesh mesh2;
  mesh2.load("square1.mesh");
  mesh2.refine_all_elements();
  mesh2.refine_all_elements();
  Solution zero;
  zero.set_const(&mesh2, 0.0);

  WeakForm wf2(1);
  wf2.add_biform(0, 0, bilinear_form2, UNSYM, 0, 1, &zero);
  wf2.add_liform(0, linear_form);

  LinSystem sys2(&wf2, &umfpack);
  sys2.set_spaces(1, &space);
  sys2.set_pss(1, &pss);


  space.set_uniform_order(5);
  space.assign_dofs();

  sys2.assemble();
  sys2.solve(1, &sln2);

  ScalarView view2("Solution 2 - laplace");
  view2.show(&sln2);
  view2.wait_for_close();

}

