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
#include "space.h"
#include "precalc.h"
#include "refmap.h"
#include "solution.h"


NonlinSystem::NonlinSystem(WeakForm* wf, Solver* solver)
            : LinSystem(wf, solver)
{
  alpha = 1.0;
  res_l2 = res_l1 = res_max = -1.0;
  // todo
  
  // tell LinSystem not to add Dirichlet contributions to the RHS
  want_dir_contrib = false;
}


void NonlinSystem::set_ic(MeshFunction* fn) 
{
  // todo
}



void NonlinSystem::assemble()
{ 
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
    
  // replace RHS with J(Y_n)*Y_n - alpha*F(Y_n)
  scalar *tmp = new scalar[ndofs];
  A_times_vec(tmp);            
  for (int i = 0; i < ndofs; i++)
    RHS[i] = tmp[i] - alpha*RHS[i]; 
  delete [] tmp;
}


void NonlinSystem::A_times_vec(scalar* result)
{
  if (!mat_sym)
  {
    if (mat_row) // CSR
    {
      for (int i = 0, n = 0; i < ndofs; i++)
        for (result[i] = 0; n < Ap[i+1]; n++)
          result[i] += Ax[n] * Vec[Ai[n]];
    }
    else // CSC
    {
      memset(result, 0, sizeof(scalar) * ndofs);
      for (int j = 0, n = 0; j < ndofs; j++)
        for ( ; n < Ap[j+1]; n++)
          result[Ai[n]] += Ax[n] * Vec[j];
    }
  }
  else
  {
    error("not implemented yet for symmetric matrices");
  }
}


