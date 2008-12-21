#include "hermes2d.h"

const double e_0   = 8.8541878176 * 1e-12;
const double mu_0   = 1.256 * 1e-6;
const double k = 1.0;

//////////////////// EXACT SOLUTION //////////////////////////////////////////////////////////////////////////

extern "C" void fresnl( double xxa, double *ssa, double *cca );


scalar Fn(double u)
{
  double s, c;
  fresnl(sqrt(2/M_PI) * u, &s , &c);
  scalar fres = complex(c,-s);
  scalar a = complex(0.0, M_PI/4);
  scalar b = complex(0.0, u*u);
  return 0.5*sqrt(M_PI) * exp(b) * (exp(-a) - sqrt(2)*(fres));
}


scalar Fder(double u)
{
  scalar a = complex(0.0, M_PI/4);
  scalar b = complex(0.0, u*u);
  scalar d = complex(0.0, 2.0*u);
  double s, c;
  fresnl(sqrt(2/M_PI) * u, &s , &c);
  scalar fres = complex(c,-s);
  scalar fresder = exp(-b);
  
  return 0.5*sqrt(M_PI) * exp(b) * ( d * (exp(-a) - sqrt(2)*(fres)) - sqrt(2)*fresder*sqrt(2.0/M_PI) );
}


scalar Fder2(double u)
{
  scalar a = complex(0.0, M_PI/4);
  scalar i = complex(0.0,1.0);
  scalar b = complex(0.0, u*u);
  scalar d = complex(0.0, 2.0*u);
  double s, c;
  fresnl(sqrt(2/M_PI) * u, &s , &c);
  scalar fres = complex(c,-s);
  scalar fresder = exp(-b);
  scalar fresder2 = exp(-b)*(-2.0 * i * u);
  
  return 2.0 * u * i * Fder(u) + 
         0.5 * sqrt(M_PI) * exp(b) * 
          ( 2.0 * i * (exp(-a) - sqrt(2)*(fres)) + d * (-sqrt(2)*fresder*sqrt(2.0/M_PI)) - sqrt(2) * fresder2 * sqrt(2.0/M_PI) );
}


scalar der_Hr(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar a = complex(0.0, M_PI/4 - k*r);
  scalar i = complex(0.0,1.0);
  return 1/sqrt(M_PI) * exp(a) * 
        ( (-i*k)*(Fn(sqrt(2*k*r)*sin(t/2 - M_PI/8)) + Fn(sqrt(2*k*r)*sin(t/2 + M_PI/8))) + 
        (Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8))*(sqrt(k)/sqrt(2*r)*sin(t/2 - M_PI/8)) + 
         Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8))*(sqrt(k)/sqrt(2*r)*sin(t/2 + M_PI/8))));
}


scalar der_Hrr(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar a = complex(0.0, M_PI/4 - k*r);
  scalar i = complex(0.0,1.0);
  scalar f1_d = Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d = Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar f1_d2 = Fder2(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d2 = Fder2(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar b1 = (sqrt(k/(2*r))*sin(t/2 - M_PI/8));
  scalar b2 = (sqrt(k/(2*r))*sin(t/2 + M_PI/8));
  return -i*k*der_Hr(x,y) + 1/sqrt(M_PI) * exp(a) * 
        ( (-i*k)*(f1_d*b1 + f2_d*b2) + 
        ( f1_d2*b1*b1 + f2_d2*b2*b2) + 
          f1_d*(-0.5*sqrt(k/(2*r*r*r))*sin(t/2 - M_PI/8))  + f2_d*(-0.5*sqrt(k/(2*r*r*r))*sin(t/2 + M_PI/8)));
}


scalar der_Hrt(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar i = complex(0.0,1.0);  
  scalar a = complex(0.0, M_PI/4 - k*r);
  scalar f1_d = Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d = Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar f1_d2 = Fder2(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d2 = Fder2(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar b1 = (sqrt(k)/sqrt(2*r)*sin(t/2 - M_PI/8));
  scalar b2 = (sqrt(k)/sqrt(2*r)*sin(t/2 + M_PI/8));
  scalar c1 = (sqrt(k*r)/sqrt(2)*cos(t/2 - M_PI/8));
  scalar c2 = (sqrt(k*r)/sqrt(2)*cos(t/2 + M_PI/8));
  return 1/sqrt(M_PI) * exp(a) * 
        ( (-i*k)*(f1_d*c1 + f2_d*c2) + 
        ( f1_d2*b1*c1 + f2_d2*b2*c2) + 
          f1_d*(0.5*sqrt(k/(2*r))*cos(t/2 - M_PI/8))  + f2_d*(0.5*sqrt(k/(2*r))*cos(t/2 + M_PI/8)));
}


scalar der_Ht(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar a = complex(0.0, M_PI/4 - k*r);  
  return 1/sqrt(M_PI) * exp(a) * 
         (Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8))*(sqrt(k*r/2)*cos(t/2 - M_PI/8)) + 
          Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8))*(sqrt(k*r/2)*cos(t/2 + M_PI/8)));
}


scalar der_Htr(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar i = complex(0.0,1.0);  
  scalar a = complex(0.0, M_PI/4 - k*r);  
  scalar f1_d = Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d = Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar f1_d2 = Fder2(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d2 = Fder2(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar b1 = (sqrt(k)/sqrt(2*r)*sin(t/2 - M_PI/8));
  scalar b2 = (sqrt(k)/sqrt(2*r)*sin(t/2 + M_PI/8));
  scalar c1 = (sqrt(k*r)/sqrt(2)*cos(t/2 - M_PI/8));
  scalar c2 = (sqrt(k*r)/sqrt(2)*cos(t/2 + M_PI/8));
  return -i*k*der_Ht(x,y) + 1/sqrt(M_PI) * exp(a) * 
         ((f1_d2*b1*c1 + f2_d2*b2*c2) + 
          f1_d*(0.5*sqrt(k/(2*r))*cos(t/2 - M_PI/8))  + f2_d*(0.5*sqrt(k/(2*r))*cos(t/2 + M_PI/8)));
}



scalar der_Htt(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar a = complex(0.0, M_PI/4 - k*r);  
  scalar f1_d = Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d = Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar f1_d2 = Fder2(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d2 = Fder2(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar c1 = (sqrt(k*r/(2))*cos(t/2 - M_PI/8));
  scalar c2 = (sqrt(k*r/(2))*cos(t/2 + M_PI/8));
  return 1/sqrt(M_PI) * exp(a) * 
         ((f1_d2*c1*c1 + f2_d2*c2*c2) + 
          f1_d*(-0.5*sqrt(k*r/2)*sin(t/2 - M_PI/8))  + f2_d*(-0.5*sqrt(k*r/2)*sin(t/2 + M_PI/8)));
}


scalar exact0(double x, double y, scalar& dx, scalar& dy)
{
  double r = sqrt(x*x + y*y);
  double theta = atan2(y,x);
  scalar Hr = der_Hr(x,y);
  scalar Ht = der_Ht(x,y);
  scalar i = complex(0.0,1.0);
  return  -i * (Hr * y/r + Ht * x/(r*r));
}


scalar exact1(double x, double y, scalar& dx, scalar& dy)
{
  double r = sqrt(x*x + y*y);
  double theta = atan2(y,x);
  scalar Hr = der_Hr(x,y);
  scalar Ht = der_Ht(x,y);
  scalar i = complex(0.0,1.0);
  return  i * ( Hr * x/r - Ht * y/(r*r));
}


static void exact_sol(double x, double y, scalar& u0, scalar& u1, scalar& u1dx, scalar& u0dy)
{
  scalar dx,dy;
  u0 = exact0(x,y,dx,dy);
  u1 = exact1(x,y,dx,dy);
  
  scalar Hr = der_Hr(x,y);
  scalar Ht = der_Ht(x,y);  
  scalar Hrr = der_Hrr(x,y);
  scalar Hrt = der_Hrt(x,y);
  scalar Htr = der_Htr(x,y);  
  scalar Htt = der_Htt(x,y);
  
  double r = sqrt(x*x + y*y);
  double theta = atan2(y,x);
  scalar i = complex(0.0,1.0);

  u1dx =  i * (( Hrr * x/r + Hrt * (-y/(r*r))) * x/r     + Hr * (y*y)/(r*r*r) - 
               ((Htr * x/r + Htt * (-y/(r*r))) * y/(r*r) + Ht * (-2.0*x*y/(r*r*r*r))));
  u0dy = -i * (( Hrr * y/r + Hrt *   x/(r*r))  * y/r     + Hr * (x*x)/(r*r*r) + 
                (Htr * y/r + Htt *   x/(r*r))  * x/(r*r) + Ht * (-2.0*x*y/(r*r*r*r)));
}


scalar2& exact(double x, double y, scalar2& dx, scalar2& dy)
{
  scalar2 ex;
  exact_sol(x,y, ex[0], ex[1], dx[1], dy[0]);

  dx[0] = 0.0; // not important
  dy[1] = 0.0; // not important

  return ex;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////// BOUNDARY CONDITIONS and BILINEAR FORM ///////////////////////////////////////////////////////

int bc_types(int marker)
{
  return BC_ESSENTIAL; 
}

double2 tau[5] = { { 0, 0}, { 1, 0 },  { 0, 1 }, { -1, 0 }, { 0, -1 } };

complex bc_values(int marker, double x, double y)
{
  scalar dx, dy;
  return exact0(x, y, dx, dy)*tau[marker][0] + exact1(x, y, dx, dy)*tau[marker][1];
}

complex bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{
  return int_curl_e_curl_f(fu, fv, ru, rv) - int_e_f(fu, fv, ru, rv);
}

////////////////////////// ADAPTIVITY /////////////////////////////////////////////////////////////////////////////

// struct ElemList
// {
//   int id;
//   double error; 
// };
// 
// 
// // Helper function for sorting elements by their error
// static int compare_error(const void* x1, const void* x2)
// {
//   const ElemList *ei1 = (const ElemList*) x1;
//   const ElemList *ei2 = (const ElemList*) x2;
//   return (ei1)->error < (ei2)->error ? 1 : -1;
// }
// 
// 
// bool adapt_solution(double eps, Mesh* mesh, Mesh* rmesh, Solution* sln, Solution* rsln)
// {
//   int i, j;
//   
//   sln->enable_transform( true);
//   
//   int ne = mesh->get_max_element_id() + 1;
//   double elem_error[ne];
//   memset(elem_error, 0, sizeof(double) * ne);
//   
//   double error;
// 
//   int num = mesh->get_num_active_elements();
//   ElemList* elist = new ElemList[num]; 
// 
//   Quad2D* quad = &g_quad_2d_std;
//   sln->set_quad_2d(quad);
//   rsln->set_quad_2d(quad);
//     
//   Mesh* meshes[2] = { mesh, rmesh };
//   Transformable* tr[2] = { sln, rsln };
//   Traverse trav;
//   trav.begin(2, meshes, tr);
// 
//   
//   Element** ee;
//   while ((ee = trav.get_next_state(NULL, NULL)) != NULL)
//   {
//     RefMap* re = sln->get_refmap();
//     RefMap* rr = rsln->get_refmap();
// 
//     error = int_hcurl_error(sln, rsln, re, rr);
//     elem_error[ee[0]->id] += error;
//   }
//   trav.finish();
//       
//   Element* e;
//   double total_err = 0.0;
//   int n = 0;
//   for_all_active_elements(e, mesh)
//   {
//     elist[n].id = e->id;
//     elist[n].error = elem_error[e->id];
//     n++;
//     
//     total_err += elem_error[e->id];
//   }
//   
//   // sort elements by their error
//   qsort(elist, n, sizeof(ElemList), compare_error);
//   
//   // refine worst elements
//   info("Refining elements:");
//   Space* space = sln->get_space();
//   for (i = 0; i < num; i++)
//   {
//     if (elist[i].error < 0.3 * elist[0].error) break;
// 
//     Element* e;
//     e = mesh->get_element(elist[i].id);
//     int p[4];
//     bool split = true;
//     int current = space->get_element_order(elist[i].id);
// 
//     rsln->set_quad_2d(&g_quad_2d_std);   
//     rsln->enable_transform(false);
//     verbose("Refining element #%d, Error %g%%", elist[i].id, elist[i].error);
//     get_optimal_refinement_curl2(e, current, rsln, split, p);
//     rsln->enable_transform(true);
//       
//     if (!split)
//       space->set_element_order(elist[i].id, p[0]);
//     else 
//     {
//       mesh->refine_element(elist[i].id);
//       for (j = 0; j < 4; j++)
//         space->set_element_order(e->sons[j]->id, p[j]);
//     }   
//   }
// 
//   delete [] elist;
//  
//   return false;
// }


int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);

  
  Mesh mesh;
  mesh.load("screen-quad.mesh");
  //mesh.load("screen-tri.mesh");

  HcurlShapesetGradLeg shapeset;
  PrecalcShapeset pss(&shapeset);

  HcurlSpace space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(1);

  DiscreteProblem ep;
  ep.set_num_equations(1);
  ep.set_spaces(1, &space);
  ep.set_pss(1, &pss);
  ep.set_bilinear_form(0, 0, NULL, bilinear_form, NULL);
  ep.set_linear_form(0, NULL, NULL);

  // Reference solution
  Mesh rmesh;
  HcurlSpace rspace(&rmesh, &shapeset);
  rspace.set_bc_types(bc_types);
  rspace.set_bc_values(bc_values);
  
  DiscreteProblem rp;
  rp.copy(&ep);
  rp.set_spaces(1, &rspace );   

  OrderView  ord("Polynomial Orders", 325, 400, 600, 600);  

  ScalarView Xview_r("Electric field X - real",  0,0,300,300);
  ScalarView Yview_r("Electric field Y - real",325,0,300,300);
  ScalarView Xview_i("Electric field X - imag",650,0,300,300);
  ScalarView Yview_i("Electric field Y - imag",975,0,300,300);
      
  GnuplotGraph graph;
  graph.set_captions("Error Convergence for the Screen Problem in H(curl)", "Degrees of Freedom", "Error [%]");
  graph.add_row("ortho adaptivity", "k", "-", "o");
  graph.set_log_y();

  double error0 = 0.0;
  int it = 1;
  begin_time();
  bool done = false;
  do
  {
    printf("\n\n---- it=%d ------------------------------------------------------------------\n\n", it++);

    space.assign_dofs();
    ep.create_matrix();
    ep.assemble_matrix_and_rhs();

    Solution sln;
    ep.solve_system(1, &sln);    
    
    ExactSolution ex(&mesh, exact);
    double error = 100 * hcurl_error(&sln, &ex);
    info("\nExact solution error: %g%%", error);
    graph.add_values(0, space.get_num_dofs(), error);
    graph.save("convergence.txt");
    
    // vizualization
    RealFilter real(&sln);
    ImagFilter imag(&sln);
    Xview_r.set_min_max_range(-3.0,1.0);
    Xview_r.show_scale(false);
    Xview_r.show(&real, 0.01 * EPS_HIGH, FN_VAL_0);
    Yview_r.set_min_max_range(-4.0,4.0);
    Yview_r.show_scale(false);
    Yview_r.show(&real, 0.01 * EPS_HIGH, FN_VAL_1);
    Xview_i.set_min_max_range(-1.0,4.0);
    Xview_i.show_scale(false);
    Xview_i.show(&imag, 0.01 * EPS_HIGH, FN_VAL_0);
    Yview_i.set_min_max_range(-4.0,4.0);
    Yview_i.show_scale(false);
    Yview_i.show(&imag, 0.01 * EPS_HIGH, FN_VAL_1);

    ord.show(&space);

    // Reference solution
    Solution rsln;
    rmesh.copy(&mesh);
    rmesh.refine_all_elements();
    rspace.copy_orders(&space, 1);
    rspace.assign_dofs();
    rp.create_matrix();
    rp.assemble_matrix_and_rhs();
    rp.solve_system(1, &rsln);
    rp.free_matrix();

    HcurlOrthoHP hp(1, &space);
    double errorestimate = hp.calc_error(&sln, &rsln) * 100;
    info("Error estimate: %g", errorestimate);
    if (errorestimate < 0.08) done = true;
    else hp.adapt(0.4, 1);

  }

  while (!done);
  verbose("\nTotal running time: %g sec", end_time());
 

  hermes2d_finalize();
  return 0;
}
