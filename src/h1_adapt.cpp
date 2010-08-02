// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "common.h"
#include "limit_order.h"
#include "forms.h"
#include "refmap.h"
#include "integrals_h1.h"
#include "element_to_refine.h"
#include "adapt.h"
#include "h1_adapt.h"
#include "linsystem.h"
#include "views/scalar_view.h"
#include "views/order_view.h"
#include "refsystem.h"

using namespace std;

H1Adapt::H1Adapt(LinSystem* ls) : Adapt(ls) 
{
  int n = ls->wf->get_neq();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == j) {
        form[i][j] = h1_form<double, scalar>;
        ord[i][j]  = h1_form<Ord, Ord>;
      }
    }
  }
}

// Mesh is adapted to represent given exact function with given accuracy
// in a given projection norm.
void adapt_to_exact_function(Space *space, ExactFunction exactfn, 
			     RefinementSelectors::Selector* selector, 
                             double threshold, int strategy, 
                             int mesh_regularity, double err_stop, 
                             int ndof_stop, int proj_norm, 
                             bool verbose, 
                             bool visualization, Solution* sln) 
{
  if (verbose == true) printf("Mesh adaptivity to given exact function:\n");

  // Initialize a dummy weak formulation.
  WeakForm wf_dummy;

  // Initialize the linear system.
  LinSystem ls(&wf_dummy, space);

  // Initialize views.
  ScalarView* view_c = new ScalarView("Coarse mesh projection", 0, 0, 410, 300);
  ScalarView* view_f = new ScalarView("Exact function on fine mesh", 420, 0, 410, 300);
  OrderView* ordview_c = new OrderView("Coarse mesh", 840, 0, 350, 300);
  OrderView* ordview_f = new OrderView("Fine mesh", 1200, 0, 350, 300);
  ScalarView* view_e = new ScalarView("Error estimate", 0, 360, 410, 300);

  // Adaptivity loop:
  Solution* sln_coarse = new Solution();
  Solution* sln_fine = new Solution();
  int as = 1; bool done = false;
  do
  {
    // Refine mesh uniformly.
    RefSystem rs(&ls);

    // Assign the function f() to the fine mesh.
    sln_fine->set_exact(rs.get_mesh(0), exactfn);

    // Project the function f() on the coarse mesh.
    ls.project_global(exactfn, sln_coarse, proj_norm);

    // Create DiffFilter for the error.
    DiffFilter difff(sln_fine, sln_coarse, H2D_FN_VAL, H2D_FN_VAL);

    // View the approximation of the exact function.
    if (visualization) {
      view_c->show(sln_coarse);
      view_f->show(sln_fine);
      char title[100];
      sprintf(title, "Coarse mesh, step %d", as);
      ordview_c->set_title(title);
      ordview_c->show(space);
      sprintf(title, "Fine mesh, step %d", as);
      ordview_f->set_title(title);
      ordview_f->show(rs.get_space(0));
      view_e->show(&difff);
      View::wait(H2DV_WAIT_KEYPRESS);
    }

    // Calculate element errors and total error estimate.
    H1Adapt hp(&ls);
    hp.set_solutions(sln_coarse, sln_fine);
    double err_est = hp.calc_error() * 100;
    if (verbose == true) printf("Step %d, ndof_coarse %d, ndof_fine %d, proj_error %g%%\n", 
                                 as, ls.get_num_dofs(), rs.get_num_dofs(), err_est);

    // If err_est too large, adapt the mesh.
    if (err_est < err_stop) done = true;
    else {
      done = hp.adapt(selector, threshold, strategy, mesh_regularity);

      if (ls.get_num_dofs() >= ndof_stop) done = true;
    }

    as++;
  }
  while (done == false);

  if (sln != NULL) sln->copy(sln_coarse);
}

// Mesh is adapted to represent a given reference solution with given accuracy
// in a given projection norm.
void adapt_to_ref_solution(Space *space, Solution* ref_sln, 
			   RefinementSelectors::Selector* selector, 
                           double threshold, int strategy, 
                           int mesh_regularity, double err_stop, 
                           int ndof_stop, int proj_norm, bool verbose, 
                           bool visualization, Solution* sln) 
{
  if (verbose == true) printf("Mesh adaptivity to given solution:\n");

  // Initialize a dummy weak formulation.
  WeakForm wf_dummy;

  // Initialize the linear system.
  LinSystem ls(&wf_dummy, space);

  // Initialize views.
  ScalarView* view_c = new ScalarView("Coarse mesh projection", 0, 0, 410, 300);
  ScalarView* view_f = new ScalarView("Solution", 420, 0, 410, 300);
  OrderView* ordview_c = new OrderView("Coarse mesh", 840, 0, 350, 300);
  OrderView* ordview_f = new OrderView("Fine mesh", 1200, 0, 350, 300);
  ScalarView* view_e = new ScalarView("Error estimate", 0, 360, 410, 300);

  // Adaptivity loop:
  Solution* sln_coarse = new Solution();
  int as = 1; bool done = false;
  do
  {
    // Refine mesh uniformly.
    RefSystem rs(&ls);

    // Project the function f() on the coarse mesh.
    ls.project_global(ref_sln, sln_coarse, proj_norm);

    // Create DiffFilter for the error.
    DiffFilter difff(ref_sln, sln_coarse, H2D_FN_VAL, H2D_FN_VAL);

    // View the approximation of the exact function.
    if (visualization) {
      view_c->show(sln_coarse);
      view_f->show(ref_sln);
      char title[100];
      sprintf(title, "Coarse mesh, step %d", as);
      ordview_c->set_title(title);
      ordview_c->show(space);
      sprintf(title, "Fine mesh, step %d", as);
      ordview_f->set_title(title);
      ordview_f->show(rs.get_space(0));
      view_e->show(&difff);
      View::wait(H2DV_WAIT_KEYPRESS);
    }

    // Calculate element errors and total error estimate.
    H1Adapt hp(&ls);
    hp.set_solutions(sln_coarse, ref_sln);
    double err_est = hp.calc_error() * 100;
    if (verbose ==  true) printf("Step %d, ndof_coarse %d, ndof_fine %d, proj_error %g%%\n", 
                                 as, ls.get_num_dofs(), rs.get_num_dofs(), err_est);

    // If err_est too large, adapt the mesh.
    if (err_est < err_stop) done = true;
    else {
      done = hp.adapt(selector, threshold, strategy, mesh_regularity);

      if (ls.get_num_dofs() >= ndof_stop) done = true;
    }

    as++;
  }
  while (done == false);

  if (sln != NULL) sln->copy(sln_coarse);
}

// Mesh is adapted to represent the Dirichlet lift with given accuracy
// in a given projection norm.
void adapt_to_dirichlet_lift(Space *space, RefinementSelectors::Selector* selector, 
                             double threshold, int strategy, 
                             int mesh_regularity, double err_stop, 
                             int ndof_stop, int proj_norm, 
                             bool use_projection, bool verbose, 
                             bool visualization, Solution* sln) 
{
  if (verbose == true) printf("Mesh adaptivity to the Dirichlet lift:\n");

  // Initialize a dummy weak formulation.
  WeakForm wf_dummy;

  // Initialize the linear system.
  LinSystem ls(&wf_dummy, space);

  // Initialize views.
  ScalarView* view_c = new ScalarView("Coarse mesh projection", 0, 0, 410, 300);
  ScalarView* view_f = new ScalarView("Solution", 420, 0, 410, 300);
  OrderView* ordview_c = new OrderView("Coarse mesh", 840, 0, 350, 300);
  OrderView* ordview_f = new OrderView("Fine mesh", 1200, 0, 350, 300);
  ScalarView* view_e = new ScalarView("Error estimate", 0, 360, 410, 300);

  // Adaptivity loop:
  Solution* sln_coarse = new Solution();
  Solution* ref_sln = new Solution();
  int as = 1; bool done = false;
  do
  {
    // Refine mesh uniformly.
    RefSystem rs(&ls);

    // Set the reference solution to be the Dirichlet lift.
    if (use_projection == true) {
      Solution dir_lift;
      dir_lift.set_dirichlet_lift(rs.get_space(0));
      rs.project_global(&dir_lift, ref_sln, proj_norm);
    }
    else ref_sln->set_dirichlet_lift(rs.get_space(0));

    // Project the reference solution on the coarse mesh.
    if (use_projection == true) {
      Solution dir_lift;
      dir_lift.set_dirichlet_lift(ls.get_space(0));
      ls.project_global(&dir_lift, sln_coarse, proj_norm);
    }
    else sln_coarse->set_dirichlet_lift(ls.get_space(0));

    // Create DiffFilter for the error.
    DiffFilter difff(ref_sln, sln_coarse, H2D_FN_VAL, H2D_FN_VAL);

    // View the approximation of the exact function.
    if (visualization) {
      view_c->show(sln_coarse);
      view_f->show(ref_sln);
      char title[100];
      sprintf(title, "Coarse mesh, step %d", as);
      ordview_c->set_title(title);
      ordview_c->show(space);
      sprintf(title, "Fine mesh, step %d", as);
      ordview_f->set_title(title);
      ordview_f->show(rs.get_space(0));
      view_e->show(&difff);
      View::wait(H2DV_WAIT_KEYPRESS);
    }

    // Calculate element errors and total error estimate.
    H1Adapt hp(&ls);
    hp.set_solutions(sln_coarse, ref_sln);
    double err_est = hp.calc_error() * 100;
    if (verbose ==  true) printf("Step %d, ndof_coarse %d, ndof_fine %d, proj_error %g%%\n", 
                                 as, ls.get_num_dofs(), rs.get_num_dofs(), err_est);

    // If err_est too large, adapt the mesh.
    if (err_est < err_stop) done = true;
    else {
      done = hp.adapt(selector, threshold, strategy, mesh_regularity);

      if (ls.get_num_dofs() >= ndof_stop) done = true;
    }

    as++;
  }
  while (done == false);

  if (sln != NULL) sln->copy(sln_coarse);
}
