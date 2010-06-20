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
#include "views/view.h"
#include "views/vector_view.h"

#include "solvers.h"

void NonlinSystem::init_nonlin()
{
  alpha = 1.0;
  //res_l2 = res_l1 = res_max = -1.0;

  // Tell LinSystem not to add Dirichlet contributions to the RHS.
  // The reason for this is that in NonlinSystem the Jacobian matrix 
  // is assembled, and the Dirichlet lift is cancelled by the derivative 
  // with respect to the coefficient vector.
  want_dir_contrib = false;
}

// this is needed because of a constructor in RefSystem
NonlinSystem::NonlinSystem() {}

NonlinSystem::NonlinSystem(WeakForm* wf_, CommonSolver* solver_)
{ 
  this->init_lin(wf_, solver_);
  this->init_nonlin();
}

NonlinSystem::NonlinSystem(WeakForm* wf_)
{
  CommonSolver *solver_ = NULL;
  this->init_lin(wf_, solver_);
  this->init_nonlin();
}

NonlinSystem::NonlinSystem(WeakForm* wf_, CommonSolver* solver_, Tuple<Space*> spaces_)
{
  int n = spaces_.size();
  if (n != wf_->neq) 
    error("Number of spaces does not match number of equations in LinSystem::LinSystem().");
  this->init_lin(wf_, solver_);
  this->init_spaces(spaces_);
  this->alloc_and_zero_vectors();
  this->init_nonlin();
}

NonlinSystem::NonlinSystem(WeakForm* wf_, Tuple<Space*> spaces_)
{
  CommonSolver* solver_ = NULL;
  this->init_lin(wf_, solver_);
  this->init_spaces(spaces_);
  this->alloc_and_zero_vectors();
  this->init_nonlin();
}

NonlinSystem::NonlinSystem(WeakForm* wf_, CommonSolver* solver_, Space *s_)
{
  if (wf_->neq != 1) 
    error("Number of spaces does not match number of equations in LinSystem::LinSystem().");
  this->init_lin(wf_, solver_);
  this->init_space(s_);
  this->alloc_and_zero_vectors();
  this->init_nonlin();
}

NonlinSystem::NonlinSystem(WeakForm* wf_, Space *s_)
{
  CommonSolver *solver_ = NULL;
  this->init_lin(wf_, solver_);
  this->init_space(s_);
  this->alloc_and_zero_vectors();
  this->init_nonlin();
}

void NonlinSystem::free()
{
  /* FIXME - MEMORY LEAK
  LinSystem::free_matrix();
  LinSystem::free_vectors();
  if (solver) solver->free_data(slv_ctx);

  struct_changed = values_changed = true;
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  wf_seq = -1;
  */
}

void NonlinSystem::assemble(bool rhsonly)
{
  // sanity checks
  int ndof = this->get_num_dofs();
  if (rhsonly) error("Parameter rhsonly = true has no meaning in NonlinSystem.");
  if (this->have_spaces == false)
    error("Before assemble(), you need to initialize spaces.");
  if (this->spaces == NULL) error("spaces = NULL in LinSystem::assemble().");

  // assemble J(Y_n) and store in A, assemble F(Y_n) and store in RHS
  LinSystem::assemble();

  // calculate norms of the residual F(Y_n)
  res_l2 = res_l1 = res_max = 0.0;
  for (int i = 0; i < ndof; i++)
  {
    res_l2 += sqr(this->RHS[i]);
    res_l1 += magn(this->RHS[i]);
    if (magn(this->RHS[i]) > res_max) res_max = magn(this->RHS[i]);
  }
  res_l2 = sqrt(res_l2);

  // multiply RHS by -alpha
  for (int i = 0; i < ndof; i++)
    this->RHS[i] *= -this->alpha;
}

// The solve() function is almost identical to the original one in LinSystem.
// It does not put the Dirichlet lift vector Dir to the right-hand side.
bool NonlinSystem::solve(Tuple<Solution*> sln)
{
  int n = sln.size();

  // if the number of solutions does not match the number of equations, throw error
  if (n != this->wf->neq) 
    error("Number of solutions does not match the number of equations in LinSystem::solve().");

  // if Vec is not initialized, throw error
  if (this->Vec == NULL) error("Vec is NULL in NonlinSystem::solve().");

  // check vector size
  int ndof = this->get_num_dofs();
  if (Vec_length != ndof || RHS_length != ndof || Dir_length != ndof)
    error("Length of vectors Vec, RHS or Dir does not match this->ndof in LinSystem::solve().");

  // check matrix size
  if (ndof == 0) error("ndof = 0 in LinSystem::solve().");
  if (ndof != this->A->get_size()) 
    error("Matrix size does not match vector length in LinSystem:solve().");

  // time measurement
  TimePeriod cpu_time;

  // solve the system - this is different from LinSystem
  scalar* delta = (scalar*) malloc(ndof * sizeof(scalar));
  memcpy(delta, this->RHS, sizeof(scalar) * ndof);
  this->solver->solve(this->A, this->Vec);
  report_time("Solved in %g s", cpu_time.tick().last());
  // add the increment dY_{n+1} to the previous solution vector
  for (int i = 0; i < ndof; i++) this->Vec[i] += delta[i];
  ::free(delta);

  // copy the solution coefficient vectors into Solutions
  cpu_time.tick(H2D_SKIP);
  for (int i = 0; i < n; i++)
  {
    sln[i]->set_fe_solution(this->spaces[i], this->pss[i], this->Vec);
  }
  report_time("Exported solution in %g s", cpu_time.tick().last());
  return true;
}

// single equation case
bool NonlinSystem::solve(Solution* sln)
{
  bool flag;
  flag = this->solve(Tuple<Solution*>(sln));
  return flag;
}

// Newton's method for an arbitrary number of equations.
bool NonlinSystem::solve_newton(Tuple<Solution*> u_prev, double newton_tol, 
                                int newton_max_iter, bool verbose, 
                                Tuple<MeshFunction*> mesh_fns) 
{
  // sanity checks
  int n = u_prev.size();
  if (n != this->wf->neq) 
    error("The number of solutions in newton_solve() must match the number of equation in the PDE system.");
  if (this->spaces == NULL) error("spaces is NULL in solve_newton().");
  for (int i=0; i < n; i++) {
    if (this->spaces[i] == NULL) error("spaces[%d] is NULL in solve_newton().", i);
  }
  int n_mesh_fns;
  if (mesh_fns == Tuple<MeshFunction*>()) n_mesh_fns = 0;
  else n_mesh_fns = mesh_fns.size();
  for (int i=0; i<n_mesh_fns; i++) {
    if (mesh_fns[i] == NULL) error("a filter is NULL in solve_newton().");
  }

  int it = 1;
  double res_l2_norm;
  do
  {
    info("---- Newton iter %d:", it); 

    // reinitialize filters
    for (int i=0; i < n_mesh_fns; i++) mesh_fns[i]->reinit();

    // assemble the Jacobian matrix and residual vector,
    // solve the system
    this->assemble();
    this->solve(u_prev);

    // calculate the l2-norm of residual vector
    res_l2_norm = this->get_residual_l2_norm();
    if (verbose) printf("---- Newton iter %d, ndof %d, res. l2 norm %g\n", 
                        it, this->get_num_dofs(), res_l2_norm);

    it++;
  }
  while (res_l2_norm > newton_tol && it <= newton_max_iter);

  // returning "true" if converged, otherwise returning "false"
  if (it <= newton_max_iter) return true;
  else return false;
}

