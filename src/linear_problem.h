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

#ifndef __H2D_LINPROBLEM_H
#define __H2D_LINPROBLEM_H

#include "matrix_old.h"
#include "filter.h"
#include "views/scalar_view.h"
#include "views/order_view.h"
#include "function.h"
#include "discrete_problem.h"

class Solution;
class MeshFunction;

H2D_API_USED_TEMPLATE(Tuple<Space*>); ///< Instantiated template. It is used to create a clean Windows DLL interface.
H2D_API_USED_TEMPLATE(Tuple<Solution*>); ///< Instantiated template. It is used to create a clean Windows DLL interface.

///
///
///
///
///
class H2D_API LinearProblem : public DiscreteProblem
{
public:
  LinearProblem();
  LinearProblem(WeakForm* _wf);
  LinearProblem(WeakForm* _wf, Space* s_);
  LinearProblem(WeakForm* wf_, Tuple<Space*> spaces_); 
  virtual ~LinearProblem();


  /// Simplified version of the assembling procedure for the user. Inside,
  /// this constructs the matrix A, vectors Dir and RHS, and adds the Dir
  /// vector to RHS.
  virtual void assemble(bool rhsonly = false);


  /// Solves the matrix problem with "mat" and "rhs", and copies the result 
  /// to the vector "vec".
  virtual bool solve(CooMatrix* mat, scalar* rhs, scalar* vec);

  /// Solves the matrix problem with this->A and this->RHS, copies the result 
  /// into this->Vec, and propagates this->Vec into one or more Solutions. 
  virtual bool solve(Tuple<Solution*> sln);
  virtual bool solve(Solution* sln);            // single equation case

  friend class RefDiscreteProblem;

};

#endif
