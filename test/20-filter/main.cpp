#include "hermes2d.h"


class TestFilter : public Filter
{
public:

  TestFilter(MeshFunction* sln1, MeshFunction* sln2)
    : Filter(sln1, sln2) {}

  virtual scalar get_pt_value(double x, double y, int item) { return 0; }

protected:
  
  virtual void precalculate(int order, int mask)
  {
    Quad2D* quad = quads[cur_quad];
    int np = quad->get_num_points(order);
    Node* node = new_node(FN_VAL_0, np);
  
    sln[0]->set_quad_order(order);
    sln[1]->set_quad_order(order);
  
    scalar* val1 = sln[0]->get_fn_values();
    scalar* val2 = sln[1]->get_fn_values();
  
    for (int i = 0; i < np; i++)
      node->values[0][0][i] = sqrt(val1[i] + val2[i]);
  
    replace_cur_node(node);
  }

};



scalar bc_values(int marker, double x, double y)
  { return x*x; }
  
scalar bilinear_form_unsym(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v(fu, fv, ru, rv); }

scalar linear_form(RealFunction* fv, RefMap* rv)
  { return -2*int_v(fv, rv); }

scalar y_square(double x, double y, scalar& dx, scalar& dy)
  { return y*y; }
  
  
  
void refine_randomly(Mesh* mesh, int steps)
{
  int i, id;
  for (i = 0; i < steps; i++)
  {
    do
      id = rand() % mesh->get_max_element_id();
    while (!mesh->get_element_fast(id)->active);
    mesh->refine_element(id, rand() % 3);
  }
}


int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);
  srand(5);
  
  Mesh mesh1;
  mesh1.load("test2.mesh");
  //refine_randomly(&mesh1, 20);
  //mesh1.refine_all_elements(1);
  mesh1.refine_towards_vertex(0, 15);
  
  H1ShapesetBeuchler shapeset;
  PrecalcShapeset pss(&shapeset);
  
  H1Space space(&mesh1, &shapeset);
  space.set_bc_values(bc_values);
  space.set_uniform_order(3);
  space.assign_dofs();

  DiscreteProblem dp;
  dp.set_num_equations(1);
  dp.set_spaces(1, &space);
  dp.set_pss(1, &pss);
  dp.set_bilinear_form(0, 0, bilinear_form_unsym);
  dp.set_linear_form(0, linear_form);

  Solution x2;
  dp.create_matrix();
  dp.assemble_matrix_and_rhs();
  dp.solve_system(1, &x2);

  ScalarView view1;
  view1.show(&x2, EPS_HIGH);

  Mesh mesh2;
  mesh2.load("test2.mesh");
  //refine_randomly(&mesh2, 20);
  //mesh2.refine_all_elements(2);
  
  Solution y2;
  y2.set_exact(&mesh2, y_square);
  ScalarView view2;
  view2.show(&y2, EPS_HIGH);
  
  TestFilter filter(&x2, &y2);
  ScalarView view3;
  view3.show(&filter, EPS_HIGH);

  hermes2d_finalize();
  return 0;
}
