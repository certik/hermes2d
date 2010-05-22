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
#include "linsystem.h"
#include "nonlinsystem.h"
#include "weakform.h"
#include "solver.h"
#include "space.h"
#include "precalc.h"
#include "refmap.h"
#include "solution.h"
#include "integrals_h1.h"

#include "python_solvers.h"

NonlinSystem::NonlinSystem(WeakForm* wf, Solver* solver)
            : LinSystem(wf, solver)
{
  alpha = 1.0;
  res_l2 = res_l1 = res_max = -1.0;

  // tell LinSystem not to add Dirichlet contributions to the RHS
  want_dir_contrib = false;
}

void NonlinSystem::free()
{
  if (RHS != NULL) { ::free(RHS); RHS = NULL; }
  if (Dir != NULL) { ::free(Dir-1); Dir = NULL; }
  if (Vec != NULL)
  {
   for (int i = 0; i < wf->neq; i++)
      if (spaces[i]->get_seq() != sp_seq[i])
        { ::free(Vec); Vec = NULL;  break; }
  }

  if (solver) solver->free_data(slv_ctx);

  struct_changed = values_changed = true;
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  wf_seq = -1;
}

void NonlinSystem::assemble(bool rhsonly)
{
  if (rhsonly) error("Parameter rhsonly = true has no meaning in NonlinSystem.");

  // assemble J(Y_n) and store in A, assemble F(Y_n) and store in RHS
  LinSystem::assemble();

  // calculate norms of the residual F(Y_n)
  res_l2 = res_l1 = res_max = 0.0;
  for (int i = 0; i < ndofs; i++)
  {
    res_l2 += sqr(RHS[i]);
    res_l1 += magn(RHS[i]);
    if (magn(RHS[i]) > res_max) res_max = magn(RHS[i]);
  }
  res_l2 = sqrt(res_l2);

  // multiply RHS by -alpha
  for (int i = 0; i < ndofs; i++)
    RHS[i] *= -alpha;
}


bool NonlinSystem::solve(int n, ...)
{
  // The solve() function is almost identical to the original one in LinSystem
  // except that Y_{n+1} = Y_{n} + dY_{n+1}
  TimePeriod cpu_time;

  // solve the system
  scalar* delta = (scalar*) malloc(ndofs * sizeof(scalar));
  memcpy(delta, RHS, sizeof(scalar) * this->A->get_size());
  solve_linear_system_scipy_umfpack(this->A, delta);
  report_time("Solved in %g s", cpu_time.tick().last());

  // if not initialized by set_ic(), assume Vec is a zero vector
  if (Vec == NULL)
  {
    Vec = (scalar*) malloc(ndofs * sizeof(scalar));
    memset(Vec, 0, ndofs * sizeof(scalar));
  }

  // add the increment dY_{n+1} to the previous solution vector
  for (int i = 0; i < ndofs; i++)
    Vec[i] += delta[i];
  ::free(delta);

  // initialize the Solution classes
  cpu_time.tick(H2D_SKIP);
  va_list ap;
  va_start(ap, n);
  if (n > wf->neq) n = wf->neq;
  for (int i = 0; i < n; i++)
  {
    Solution* sln = va_arg(ap, Solution*);
    sln->set_fe_solution(spaces[i], pss[i], Vec);
  }
  va_end(ap);
  report_time("Exported solution in %g s", cpu_time.tick().last());

  return true;
}

// Newton's loop for one equation
bool NonlinSystem::solve_newton(Solution* u_prev, double newton_tol, int newton_max_iter,
                                  Filter* f1, Filter* f2, Filter* f3) {
    int it = 1;
    double res_l2_norm;
    Solution sln_iter;
    Space *space = this->get_space(0);
    do
    {
      info("---- Newton iter %d:", it); it++;
      //printf("ndof = %d\n", space->get_num_dofs());

      // reinitialization of filters (if relevant)
      if (f1 != NULL) f1->reinit();
      if (f2 != NULL) f2->reinit();
      if (f3 != NULL) f3->reinit();

      // assemble the Jacobian matrix and residual vector,
      // solve the system
      this->assemble();
      this->solve(1, &sln_iter);

      // calculate the l2-norm of residual vector
      res_l2_norm = this->get_residual_l2_norm();
      info("Residual L2 norm: %g", res_l2_norm);

      // save the new solution as "previous" for the
      // next Newton's iteration
      u_prev->copy(&sln_iter);
    }
    while (res_l2_norm > newton_tol && it <= newton_max_iter);

    // returning "true" if converged, otherwise returning "false"
    if (it <= newton_max_iter) return true;
    else return false;
}

// Newton's loop for two equations
bool NonlinSystem::solve_newton(Solution* u_prev_1, Solution* u_prev_2, double newton_tol, int newton_max_iter,
                                  Filter* f1, Filter* f2, Filter* f3) {
    int it = 1;
    double res_l2_norm;
    Solution sln_iter_1, sln_iter_2;
    Space *space_1 = this->get_space(0);
    Space *space_2 = this->get_space(1);
    do
    {
      info("---- Newton iter %d:", it); it++;
      //printf("ndof = %d\n", space_1->get_num_dofs() + space_2->get_num_dofs());

      // reinitialization of filters (if relevant)
      if (f1 != NULL) f1->reinit();
      if (f2 != NULL) f2->reinit();
      if (f3 != NULL) f3->reinit();

      // assemble the Jacobian matrix and residual vector,
      // solve the system
      this->assemble();
      this->solve(2, &sln_iter_1, &sln_iter_2);

      // calculate the l2-norm of residual vector
      res_l2_norm = this->get_residual_l2_norm();
      info("Residual L2 norm: %g", res_l2_norm);

      // save the new solutions as "previous" for the
      // next Newton's iteration
      u_prev_1->copy(&sln_iter_1);
      u_prev_2->copy(&sln_iter_2);
    }
    while (res_l2_norm > newton_tol && it <= newton_max_iter);

    // returning "true" if converged, otherwise returning "false"
    if (it <= newton_max_iter) return true;
    else return false;
}

// Newton's loop for three equations
bool NonlinSystem::solve_newton(Solution* u_prev_1, Solution* u_prev_2, Solution* u_prev_3,
                                  double newton_tol, int newton_max_iter,
                                  Filter* f1, Filter* f2, Filter* f3) {
    int it = 1;
    double res_l2_norm;
    Solution sln_iter_1, sln_iter_2, sln_iter_3;
    Space *space_1 = this->get_space(0);
    Space *space_2 = this->get_space(1);
    Space *space_3 = this->get_space(2);
    do
    {
      info("---- Newton iter %d:", it); it++;
      //printf("ndof = %d\n", space_1->get_num_dofs() + space_2->get_num_dofs() + space_3->get_num_dofs());

      // reinitialization of filters (if relevant)
      if (f1 != NULL) f1->reinit();
      if (f2 != NULL) f2->reinit();
      if (f3 != NULL) f3->reinit();

      // assemble the Jacobian matrix and residual vector,
      // solve the system
      this->assemble();
      this->solve(3, &sln_iter_1, &sln_iter_2, &sln_iter_3);

      // calculate the l2-norm of residual vector
      res_l2_norm = this->get_residual_l2_norm();
      info("Residual L2 norm: %g", res_l2_norm);

      // save the new solutions as "previous" for the
      // next Newton's iteration
      u_prev_1->copy(&sln_iter_1);
      u_prev_2->copy(&sln_iter_2);
      u_prev_3->copy(&sln_iter_3);
    }
    while (res_l2_norm > newton_tol && it <= newton_max_iter);

    // returning "true" if converged, otherwise returning "false"
    if (it <= newton_max_iter) return true;
    else return false;
}
