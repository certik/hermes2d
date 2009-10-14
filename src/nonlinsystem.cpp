// This file is part of Hermes2D.
//
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@unr.edu>
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
#include "linsystem.h"
#include "nonlinsystem.h"
#include "weakform.h"
#include "solver.h"
#include "space.h"
#include "precalc.h"
#include "refmap.h"
#include "solution.h"
#include "integrals_h1.h"


NonlinSystem::NonlinSystem(WeakForm* wf, Solver* solver)
            : LinSystem(wf, solver)
{
  alpha = 1.0;
  res_l2 = res_l1 = res_max = -1.0;

  // tell LinSystem not to add Dirichlet contributions to the RHS
  want_dir_contrib = false;
}


static MeshFunction* tmp_w;

static scalar H1projection_biform(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{  return int_grad_u_grad_v(fu, fv, ru, rv) + int_u_v(fu, fv, ru, rv); }

static scalar H1projection_liform(RealFunction* fv, RefMap* rv)
{  return int_grad_w_grad_v(tmp_w, fv, rv) + int_w_v(tmp_w, fv, rv); }

static scalar L2projection_biform(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
{  return int_u_v(fu, fv, ru, rv); }

static scalar L2projection_liform(RealFunction* fv, RefMap* rv)
{  return int_w_v(tmp_w, fv, rv); }

void NonlinSystem::set_ic(MeshFunction* fn, Solution* result, int proj_norm)
{
  if (!have_spaces)
    error("You need to call set_ic() after calling set_spaces().");

  WeakForm* wf_orig = wf;

  WeakForm wf_proj(1);
  wf = &wf_proj;
  if (proj_norm == 1)  {
    wf->add_biform(0, 0, H1projection_biform);
    wf->add_liform(0, H1projection_liform, ANY, 1, fn);
  }
  else  {
    wf->add_biform(0, 0, L2projection_biform);
    wf->add_liform(0, L2projection_liform, ANY, 1, fn);
  }

  want_dir_contrib = true;
  tmp_w = fn;
  LinSystem::assemble();
  LinSystem::solve(1, result);
  want_dir_contrib = false;

  wf = wf_orig;

}

void NonlinSystem::set_vec_zero()
{
  if (Vec != NULL) ::free(Vec);
  Vec = (scalar*) malloc(ndofs * sizeof(scalar));
  for(int i=0; i<ndofs; i++) Vec[i] = 0;
}




void NonlinSystem::assemble(bool rhsonly)
{
  if (rhsonly) error("rhsonly has no meaning in NonlinSystem.");

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
  begin_time();

  // perform symbolic analysis of the matrix
  if (struct_changed)
  {
    solver->analyze(slv_ctx, ndofs, Ap, Ai, Ax, false);
    struct_changed = false;
  }

  // factorize the stiffness matrix, if needed
  if (struct_changed || values_changed)
  {
    solver->factorize(slv_ctx, ndofs, Ap, Ai, Ax, false);
    values_changed = false;
  }

  // solve the system
  scalar* delta = (scalar*) malloc(ndofs * sizeof(scalar));
  solver->solve(slv_ctx, ndofs, Ap, Ai, Ax, false, RHS, delta);
  verbose("  (total solve time: %g sec)", end_time());

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
  begin_time();
  va_list ap;
  va_start(ap, n);
  if (n > wf->neq) n = wf->neq;
  for (int i = 0; i < n; i++)
  {
    Solution* sln = va_arg(ap, Solution*);
    sln->set_fe_solution(spaces[i], pss[i], Vec);
  }
  va_end(ap);
  verbose("Exported solution in %g sec", end_time());

  return true;
}
