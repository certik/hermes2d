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

// Mesh is adapted to represent initial condition with given accuracy
// in a given projection norm.
void adapt_to_exact_function(Space *space, ExactFunction exactfn, 
				RefinementSelectors::Selector* selector, 
                                double threshold, int strategy, 
                                int mesh_regularity, double err_stop, 
                                int ndof_stop, int proj_norm, 
                                bool project_on_fine_mesh, bool verbose, 
                                Solution* sln) 
{
  if (verbose == true) printf("Mesh adaptivity to an exact function:\n");

  // Initialize a dummy weak formulation.
  WeakForm wf_dummy;

  // Initialize the linear system.
  LinSystem ls(&wf_dummy, space);
  if (verbose) printf("ndof_coarse = %d\n", ls.get_num_dofs());

  // Initialize views.
  ScalarView* view_c = new ScalarView("Coarse mesh projection", 0, 0, 410, 300);
  ScalarView* view_f = new ScalarView("Fine mesh projection", 420, 0, 410, 300);
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
    if (verbose) printf("ndof_fine = %d\n", rs.get_num_dofs());

    // Assign the function f() to the fine mesh exactly or use global projection.
    if (project_on_fine_mesh == true) rs.project_global(exactfn, sln_fine, 1);
    else sln_fine->set_exact(rs.get_mesh(0), exactfn);

    // Project the function f() on the coarse mesh.
    ls.project_global(exactfn, sln_coarse, 1);

    // Create DiffFilter for the error.
    DiffFilter difff(sln_fine, sln_coarse, H2D_FN_VAL, H2D_FN_VAL);

    // View the approximation of the exact function.
    if (verbose) {
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
      //View::wait(H2DV_WAIT_KEYPRESS);
    }

    // Calculate element errors and total error estimate.
    H1Adapt hp(&ls);
    hp.set_solutions(sln_coarse, sln_fine);
    double err_est = hp.calc_error() * 100;
    if (verbose ==  true) printf("Step %d, ndof %d, proj_error %g%%\n", 
                                 as, ls.get_num_dofs(), err_est);

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
