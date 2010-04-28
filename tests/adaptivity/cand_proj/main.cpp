#include <sstream>
#include <list>
#include <hermes2d.h>
#include <solver_umfpack.h>
#include "functions.h"

// This test tests projection of a candidate in H1 space on a quad.

/* global definitions */
#undef ERROR_SUCCESS
#undef ERROR_FAILURE
#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
#define H2D_TEST_P_INIT 1 /* initial order on a mesh */
#define H2D_TEST_ELEM_ID 0 /* ID of an alement which is going to be handled */
#define H2D_TEST_ZERO 1e-13 /* Numerical zero for double datatype */
#define delete_not_null(__ptr) if (__ptr != NULL) delete __ptr;

/* global structures */
struct TestCase {
  std::string title;
  ValueFunction func_val, func_dx, func_dy;
  int func_quad_order, start_quad_order;
  int res_split;
  int4 res_orders;

  TestCase(std::string title
    , ValueFunction func_val, ValueFunction func_dx, ValueFunction func_dy
    , int start_quad_order, int func_quad_order
    , int res_split, int res_order0, int res_order1 = 0, int res_order2 = 0, int res_order3 = 0)
    : title(title), func_val(func_val), func_dx(func_dx), func_dy(func_dy)
    , start_quad_order(start_quad_order), func_quad_order(func_quad_order), res_split(res_split) {
      res_orders[0] = res_order0;
      res_orders[1] = res_order1;
      res_orders[2] = res_order2;
      res_orders[3] = res_order3;
  };
  bool is_match(const ElementToRefine& refin) { ///< Returns true in a case of a match.
    if (res_split != refin.split) {
      info("split type %d does not match, expected %d", res_split, refin.split);
      return false;
    }
    int num_sons = refin.get_num_sons();
    for(int i = 0; i < num_sons; i++) {
      if (res_orders[i] != refin.p[i]) {
        info("order (%d, %d) of son %d does not match, expected (%d, %d)"
          , H2D_GET_H_ORDER(refin.p[i]), H2D_GET_V_ORDER(refin.p[i])
          , i
          , H2D_GET_H_ORDER(res_orders[i]), H2D_GET_V_ORDER(res_orders[i]));
        return false;
      }
    }
    return true;
  };
};

/* global variables */
Mesh* mesh = NULL;
H1Shapeset* shapeset = NULL;
PrecalcShapeset* pss = NULL;
H1Space* space = NULL;
TestCase* cur_test_case = NULL;
WeakForm* weakform = NULL;
UmfpackSolver* solver = NULL;

std::list<TestCase> test_cases;

/// Boundary condition types
int bc_types(int marker) { return BC_NONE; }

/// Dirichlet boundary condition values
scalar bc_values(int marker, double x, double y) { return 0.0; }

/// Bilinear form.
template<typename T, typename D>
D biform(int n, double *wt, Func<T> *u, Func<T> *v, Geom<T> *e, ExtData<D> *data) {
  return int_u_v<T, T>(n, wt, u, v) + int_grad_u_grad_v<T, T>(n, wt, u, v);
}

/// Linear form: order
Ord liform(int point_cnt, double *weights, Func<Ord> *values_v, Geom<Ord> *geometry, ExtData<Ord> *values_fnc_ext) {
  return Ord(H2D_GET_H_ORDER(cur_test_case->func_quad_order) + H2D_GET_V_ORDER(cur_test_case->func_quad_order) + 2*values_v->val->get_order());
}

/// Linear form: value
scalar liform(int point_cnt, double *weights, Func<double> *values_v, Geom<double> *geometry, ExtData<scalar> *values_fnc_ext) {
  scalar result = 0;
  for(int i = 0; i < point_cnt; i++) {
    double x = geometry->x[i], y = geometry->y[i];
    result += weights[i] * (values_v->val[i] * cur_test_case->func_val(x, y)
      + values_v->dx[i] * cur_test_case->func_dx(x, y)
      + values_v->dy[i] * cur_test_case->func_dy(x, y));
  }
  return result;
}

/// Initialize internal data structures.
bool init(bool tri) {
  try {
    // shapeset and cache
    shapeset = new H1Shapeset();
    pss = new PrecalcShapeset(shapeset);

    // mesh (1 element)
    mesh = new Mesh();
    if (tri) {
      const int vertex_num = 3, tria_num = 1, quad_num = 0, marker_num = 3;
      double2 vertex_array[3] = {{-1,-1}, { 1,-1}, {-1, 1}};
      int4 tria_array[1] = {{0, 1, 2, 0}};
      int5 *quad_array = NULL;
      int3 marker_array[3] = {{0, 1, 2}, {1, 2, 1}, {2, 0, 3}};
      mesh->create(vertex_num, vertex_array, tria_num, tria_array, quad_num, quad_array, marker_num, marker_array);
    }
    else {
      const int vertex_num = 4, tria_num = 0, quad_num = 1, marker_num = 4;
      double2 vertex_array[4] = {{-1,-1}, { 1,-1}, { 1, 1}, {-1, 1}};
      int4 *tria_array = NULL;
      int5 quad_array[1] = {{0, 1, 2, 3, 0}};
      int3 marker_array[4] = {{0, 1, 3}, {1, 2, 2}, {2, 3, 4}, {3, 0, 1}};
      mesh->create(vertex_num, vertex_array, tria_num, tria_array, quad_num, quad_array, marker_num, marker_array);
    }

    // finite element space
    space = new H1Space(mesh, shapeset);
    space->set_bc_types(bc_types);
    space->set_bc_values(bc_values);
    space->set_uniform_order(H2D_TEST_P_INIT);
    space->assign_dofs();

    // weakform
    weakform = new WeakForm(1);
    weakform->add_biform(0, 0, callback(biform), H2D_SYM, 0);
    weakform->add_liform(0, liform, liform, H2D_ANY, 0);

    //solver
    solver = new UmfpackSolver();

    //test cases
    if (tri) {
      test_cases.push_back(TestCase("x^2 y^2", func_x2y2_dx, func_x2y2_dx, func_x2y2_dy, 1, 2+2, H2D_REFINEMENT_P, H2D_MAKE_QUAD_ORDER(2, 0)));
    }
    else {
      test_cases.push_back(TestCase("x^2", func_x2_val, func_x2_dx, func_x2_dy, H2D_MAKE_QUAD_ORDER(1,1), H2D_MAKE_QUAD_ORDER(2,0), H2D_REFINEMENT_P, H2D_MAKE_QUAD_ORDER(2,1)));
      //test_cases.push_back(TestCase("abs(y)", func_absy_val, func_absy_dx, func_absy_dy, H2D_MAKE_QUAD_ORDER(1,1), H2D_MAKE_QUAD_ORDER(0,1), H2D_REFINEMENT_ANISO_H, H2D_MAKE_QUAD_ORDER(1,1), H2D_MAKE_QUAD_ORDER(1,1)));
      //test_cases.push_back(TestCase("abs(x)*abs(y)", func_absx_absy_val, func_absx_absy_dx, func_absx_absy_dy, H2D_MAKE_QUAD_ORDER(1,1), H2D_MAKE_QUAD_ORDER(1,1), H2D_REFINEMENT_H, H2D_MAKE_QUAD_ORDER(1,1), H2D_MAKE_QUAD_ORDER(1,1), H2D_MAKE_QUAD_ORDER(1, 1), H2D_MAKE_QUAD_ORDER(1,1)));
      test_cases.push_back(TestCase("x^2 y^2", func_x2y2_val, func_x2y2_dx, func_x2y2_dy, H2D_MAKE_QUAD_ORDER(1,1), H2D_MAKE_QUAD_ORDER(2,2), H2D_REFINEMENT_P, H2D_MAKE_QUAD_ORDER(2,2)));
      test_cases.push_back(TestCase("x^3 y", func_x3y1_val, func_x3y1_dx, func_x3y1_dy, H2D_MAKE_QUAD_ORDER(2,2), H2D_MAKE_QUAD_ORDER(3,1), H2D_REFINEMENT_P, H2D_MAKE_QUAD_ORDER(3,2)));
      test_cases.push_back(TestCase("x^3 y^4", func_x3y4_val, func_x3y4_dx, func_x3y4_dy, H2D_MAKE_QUAD_ORDER(3,3), H2D_MAKE_QUAD_ORDER(3,4), H2D_REFINEMENT_P, H2D_MAKE_QUAD_ORDER(3,4)));
    }

    return true;
  }
  catch (std::exception& e) {
    info("failed: %s", e.what());
    return false;
  }
}

/// Cleanup
void cleanup() {
  delete_not_null(solver);
  delete_not_null(weakform);
  delete_not_null(space);
  delete_not_null(mesh);
  delete_not_null(pss);
  delete_not_null(shapeset);
}

/// Test
int test() {
   bool failed = false;

  //prepare selector
  RefinementSelectors::H1NonUniformHP selector(false, RefinementSelectors::H2DRS_CAND_HP, 1.0, H2DRS_DEFAULT_ORDER, shapeset);
  std::list<TestCase>::iterator iter = test_cases.begin();
  while(iter != test_cases.end()){
    //set a current function and print info
    cur_test_case = &(*iter);
    info("test case: %s", cur_test_case->title.c_str());

    //change mesh
    space->set_element_order(H2D_TEST_ELEM_ID, cur_test_case->start_quad_order);

    //create and solve the reference system
    Solution sln, rsln;
    LinSystem ls(weakform, solver);
    ls.set_spaces(1, space);
    ls.set_pss(1, pss);
    ls.assemble();
    ls.solve(1, &sln);
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &rsln);

    //select candidate
    ElementToRefine refinement;
    Element* e = mesh->get_element(H2D_TEST_ELEM_ID);
    int order = space->get_element_order(H2D_TEST_ELEM_ID);
    selector.select_refinement(e, order, &rsln, refinement);

    //check selected candidate
    if (cur_test_case->is_match(refinement)) {
      info("  selected candidate: correct");

      //check if all candidates with higher orders has zero error
      int min_order_h = H2D_GET_H_ORDER(cur_test_case->func_quad_order);
      int min_order_v = H2D_GET_V_ORDER(cur_test_case->func_quad_order);
      const std::vector<RefinementSelectors::H1UniformHP::Cand>& candidates = selector.get_candidates();
      std::vector<RefinementSelectors::H1UniformHP::Cand>::const_iterator cand = candidates.begin();
      while (cand != candidates.end()) {
        //check if candidate qualifies
        int num_sons = cand->get_num_sons();
        bool valid = true;
        for(int i = 0; i < num_sons; i++) {
          if (H2D_GET_H_ORDER(cand->p[i]) < min_order_h || H2D_GET_V_ORDER(cand->p[i]) < min_order_v)
            valid = false;
        }

        //check if it fails or not
        if (valid && cand->error > H2D_TEST_ZERO) {
          std::stringstream str;
          str << "  " << *cand;
          info("%s", str.str().c_str());
          failed = true;
        }

        //next candidate
        cand++;
      }
      if (!failed) {
        info("  candidate projection: correct");
      }
    }
    else
      failed = true;

    //next test case
    iter++;
  }

  if (failed)
    return ERROR_FAILURE;
  else
    return ERROR_SUCCESS;
}


int main(int argc, char* argv[]) {
  //check if triangle is requested
  bool tri = false;
  if (argc > 1 && strcmp(argv[1], "-tri") == 0)
    tri = true;

  //initialize
  int result = ERROR_FAILURE;
  if (init(tri)) {
    result = test();
  }
  cleanup();
  return result;
}
