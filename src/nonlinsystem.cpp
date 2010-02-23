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


template<typename Real, typename Scalar>
Scalar H1projection_biform(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar H1projection_liform(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i] + ext->fn[0]->dx[i] * v->dx[i] + ext->fn[0]->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar L2projection_biform(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar L2projection_liform(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i]);
  return result;
}


void NonlinSystem::free()
{
  if (Ap != NULL) { ::free(Ap); Ap = NULL; }
  if (Ai != NULL) { ::free(Ai); Ai = NULL; }
  if (Ax != NULL) { ::free(Ax); Ax = NULL; }
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


void NonlinSystem::set_ic_n(int proj_norm, int n, ...)
{
  if (!have_spaces)
    error("You need to call set_ic() after calling set_spaces().");
  if (n != wf->neq || n > 10)
    error("Wrong number of initial conditions.");

  int i;
  MeshFunction* fn[10];
  Solution* result[10];

  va_list ap;
  va_start(ap, n);
  for (i = 0; i < n; i++)
    fn[i] = va_arg(ap, MeshFunction*);
  for (i = 0; i < n; i++)
    result[i] = va_arg(ap, Solution*);
  va_end(ap);

  WeakForm* wf_orig = wf;

  WeakForm wf_proj(n);
  wf = &wf_proj;
  if (proj_norm == 1)
    for (i = 0; i < n; i++)
    {
      wf->add_biform(i, i, H1projection_biform<double, scalar>, H1projection_biform<Ord, Ord>);
      wf->add_liform(i, H1projection_liform<double, scalar>, H1projection_liform<Ord, Ord>, ANY, 1, fn[i]);
    }
  else
    for (i = 0; i < n; i++)
    {
      wf->add_biform(i, i, L2projection_biform<double, scalar>, L2projection_biform<Ord, Ord>);
      wf->add_liform(i, L2projection_liform<double, scalar>, L2projection_liform<Ord, Ord>, ANY, 1, fn[i]);
    }

  want_dir_contrib = true;
  LinSystem::assemble();
  LinSystem::solve(n, result[0], result[1], result[2], result[3], result[4],
                      result[5], result[6], result[7], result[8], result[9]);
  want_dir_contrib = false;

  wf = wf_orig;
  wf_seq = -1;
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
