// This file is part of Hermes2D.
//
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@utep.edu>
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
#include "traverse.h"
#include "space.h"
#include "precalc.h"
#include "refmap.h"
#include "solution.h"

NonlinSystem::NonlinSystem(WeakForm* wf, Solver* solver) : LinSystem(wf, solver)
{
  Y = NULL;
  newton_ic_defined = false;
}

void NonlinSystem::set_newton_ic_zero(int n) 
{   
  if (Y != NULL) delete [] Y;
  Y = new scalar[n];
  for (int i = 0; i < n; i++) Y[i] = 0;
  newton_ic_defined = true; 
}

void NonlinSystem::set_newton_ic_const(int n, scalar const_val) 
{ 
  if (Y != NULL) delete [] Y;
  Y = new scalar[n];
  for (int i = 0; i < n; i++) Y[i] = const_val;
  newton_ic_defined = true; 
}

void NonlinSystem::set_newton_ic_function(int n, RealFunction *f) 
{ error("set_newton_ic_function() not implemented yet."); 
  if (Y != NULL) delete [] Y;
  Y = new scalar[n];
  scalar to_be_defined;
  for (int i = 0; i < n; i++) Y[i] = to_be_defined;
  newton_ic_defined = true; 
}

void NonlinSystem::A_times_Y(scalar *vec) 
{
  if(mat_sym) warn("Economical handling of symmetric matrices not implemented yet.");

  // TEMPORARY VERSION
  if (mat_row) A_times_Y_csr(vec);
  else A_times_Y_csc(vec);

  /* CORRECT VERSION - TO BE USED IN THE FUTURE
  if (mat_sym) {
    if (mat_row) A_times_Y_csr_sym(vec);
    else A_times_Y_csc_sym(vec);
  }
  else {
    if (mat_row) A_times_Y_csr(vec);
    else A_times_Y_csc(vec);
  }
  */
}

void NonlinSystem::A_times_Y_csr(scalar *vec) 
{ int count = 0;
  for (int row = 0; row < ndofs; row++) {
    vec[row] = 0;
    int nnz = Ap[row+1] - Ap[row]; 
    for (int c = 0; c < nnz; c++) {   
      vec[row] += Ax[count] * Y[Ai[count]];
      count++;
    }
  }
}

void NonlinSystem::A_times_Y_csc(scalar *vec) 
{ int count = 0;
  for (int row = 0; row < ndofs; row++) vec[row] = 0;
  for (int col = 0; col < ndofs; col++) {
    int nnz = Ap[col+1] - Ap[col]; 
    for (int r = 0; r < nnz; r++) {   
      vec[Ai[count]] += Ax[count] * Y[col];
      count++;
    }
  }
}

void NonlinSystem::A_times_Y_csr_sym(scalar *vec) 
{ error("A_times_Y_csr_sym() not implemented yet."); }

void NonlinSystem::A_times_Y_csc_sym(scalar *vec) 
{ error("A_times_Y_csc_sym() not implemented yet."); }

void NonlinSystem::assemble_newton(double alpha) 
{ 
  // sanity check
  if (!newton_ic_defined || Y == NULL) 
    error("Initial condition for Newton's iteration not defined.");

  // assemble J(Y_n) and store in A, assemble F(Y_n) and store in RHS
  assemble(false);

  // delete Dirichlet contributions made in
  // LinSystem:assemble() 
  for (int i = 0; i < ndofs; i++) {
    RHS[i] -= Dir[i];  
  }

  //replace RHS with J(Y_n) times Y_n minus alpha*F(Y_n)
  scalar *vec = new scalar[ndofs];
  A_times_Y(vec);            
  for (int i=0; i<ndofs; i++) RHS[i] = vec[i] - alpha*RHS[i]; 
  delete [] vec;
}

bool NonlinSystem::solve_newton(int n, ...) 
{ // This is almost identical to LinSystem::Solve(). The only exception is 
  // that the coefficient vector is saved to Y rather than to a temporary vector "vec"
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
  // scalar* vec = new scalar[ndofs];
  solver->solve(slv_ctx, ndofs, Ap, Ai, Ax, false, RHS, Y);
  verbose("  (total solve time: %g sec)", end_time());
      
  // initialize the Solution classes
  begin_time();
  va_list ap;
  va_start(ap, n);
  if (n > wf->neq) n = wf->neq;
  for (int i = 0; i < n; i++)
  {
    Solution* sln = va_arg(ap, Solution*);
    sln->set_fe_solution(spaces[i], pss[i], Y);
  }
  va_end(ap);
  verbose("Exported solution in %g sec", end_time());
  
  // delete [] vec;
  return true;
}

void NonlinSystem::compute_residuum_vector(scalar*& res) // calculates the residuum F(Y_{n+1})
{ assemble(true); 
  // delete Dirichlet contributions made in
  // LinSystem:assemble() 
  for (int i = 0; i < ndofs; i++) {
    RHS[i] -= Dir[i];  
  }
  res = this->RHS; 
}

double NonlinSystem::compute_residuum_l2_norm() 
{ scalar *res; 
  compute_residuum_vector(res); 
  double norm = 0;
  for (int i=0; i<ndofs; i++) norm += sqr(res[i]);
  return sqrt(norm); 
}

double NonlinSystem::compute_residuum_l1_norm() 
{ scalar *res; 
  compute_residuum_vector(res); 
  double norm = 0;
  for (int i=0; i<ndofs; i++) norm += magn(res[i]);
  return norm; 
}

double NonlinSystem::compute_residuum_max_norm() 
{ scalar *res; 
  compute_residuum_vector(res); 
  double norm = 0;
  for (int i=0; i<ndofs; i++) if (magn(res[i]) > norm) norm = magn(res[i]);
  return norm; 
}

void NonlinSystem::free()
{
  if (Y != NULL) { ::free(Y); Y = NULL; }
  LinSystem::free();
}
