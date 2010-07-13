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
void adapt_to_exact_function_h1(Space *space, ExactFunction exactfn, 
				RefinementSelectors::Selector* selector, double threshold, int strategy, 
                                int mesh_regularity, double err_stop, int ndof_stop, bool verbose,
                                Solution* sln) 
{
  if (verbose == true) printf("Mesh adaptivity to an exact function:\n");

  // Initialize a dummy weak formulation.
  WeakForm wf_dummy;

  // Initialize the linear system.
  LinSystem ls(&wf_dummy, space);

  // Initialize views.
  ScalarView* view = new ScalarView("Projection of initial condition", 0, 0, 410, 300);
  OrderView* ordview = new OrderView("Initial mesh", 420, 0, 350, 300);
  view->fix_scale_width(80);

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
    ls.project_global(exactfn, sln_coarse);

    // Calculate element errors and total error estimate.
    H1Adapt hp(&ls);
    hp.set_solutions(sln_coarse, sln_fine);
    double err_est = hp.calc_error() * 100;
    if (verbose ==  true) printf("Step %d, ndof %d, proj_error %g%%\n", 
                 as, ls.get_num_dofs(), err_est);

    // If err_est too large, adapt the mesh.
    if (err_est < err_stop) done = true;
    else {
      double to_be_processed = 0;
      done = hp.adapt(selector, threshold, strategy, mesh_regularity, to_be_processed);

      if (ls.get_num_dofs() >= ndof_stop) done = true;
    }

    // View the approximation of the exact function.
    if (verbose) {
      view->show(sln_coarse);
      char title[100];
      sprintf(title, "Initial mesh, step %d", as);
      ordview->set_title(title);
      ordview->show(space);
      //View::wait(H2DV_WAIT_KEYPRESS);
    }

    as++;
  }
  while (done == false);

  if (sln != NULL) sln->copy(sln_coarse);
}
