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

#include "common.h"
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

  /// Version for linear problems -- adds the dir vector to rhs.
  virtual void assemble(Matrix* mat_ext, Vector* rhs_ext, bool rhsonly = false);

  friend class RefDiscreteProblem;

};

void init_matrix_solver(MatrixSolverType matrix_solver, int ndof, 
                        Matrix* &mat, Vector* &vec, CommonSolver* &solver);

bool solve_linear(Tuple<Space *> spaces, WeakForm* wf, Tuple<Solution *> solutions, 
                  MatrixSolverType matrix_solver);

#endif
