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

// $Id: linsystem.h 1037 2008-10-01 21:32:06Z jakub $

#ifndef __HERMES2D_NONLINSYSTEM_H
#define __HERMES2D_NONLINSYSTEM_H

#include "matrix.h"
#include "function.h"

class Solution;

///  
///
///
///
///
class NonlinSystem : public LinSystem 
{
public:
 
  NonlinSystem(WeakForm* wf, Solver* solver);

  void set_newton_ic_zero(int n);  // initial condition must be defined before the Newton's iteration is run
  void set_newton_ic_const(int n, scalar const_val);                   
  void set_newton_ic_function(int n, RealFunction *f); 
  void A_times_Y(scalar *vec);     // multiply A with Y and save result in vec
  void assemble_newton(double alpha);          // assemble J(Y_n) and store in A, 
                                               // assemble F(Y_n) and store 
                                               // J(Y_n)Y_n - alpha*F(Y_n) in RHS
  bool solve_newton(int n, ...);               // solves the linear system and saves result in Y 
  void compute_residuum_vector(scalar*& res);  // calculates the residuum F(Y_{n+1})
  double compute_residuum_l2_norm();           // calculates the residuum's l2-norm 
  double compute_residuum_l1_norm();           // calculates the residuum's l1-norm 
  double compute_residuum_max_norm();          // calculates the residuum's max-norm 
  virtual void free();

protected:

  bool newton_ic_defined; 

  scalar* Y;                     ///< vector of solution coefficients (Y_n)

  void A_times_Y_csr(scalar *vec);
  void A_times_Y_csr_sym(scalar *vec);
  void A_times_Y_csc(scalar *vec);
  void A_times_Y_csc_sym(scalar *vec);
};


#endif
