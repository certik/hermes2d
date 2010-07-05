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
#include "discrete_problem.h"
#include "linear_problem.h"
#include "weakform.h"
#include "solver.h"
#include "space.h"
#include "precalc.h"
#include "refmap.h"
#include "solution.h"
#include "integrals_h1.h"
#include "views/view.h"
#include "views/vector_view.h"

LinearProblem::LinearProblem() : DiscreteProblem() {};
LinearProblem::LinearProblem(WeakForm* wf_, Space* s_) : DiscreteProblem(wf_, s_) {};
LinearProblem::LinearProblem(WeakForm* wf_, Tuple<Space*> spaces_) : DiscreteProblem(wf_, spaces_) {};
LinearProblem::~LinearProblem() {};

void LinearProblem::assemble(bool rhsonly)
{
  // sanity checks
  int ndof = this->get_num_dofs();

  if (this->have_spaces == false)
    error("Before assemble(), you need to initialize spaces.");
  if (this->spaces == NULL) error("spaces = NULL in LinearProblem::assemble().");

  // Assemble the matrix A and vector RHS. If the problem is linear,
  // we need to subtract the vector Dir from RHS.
  DiscreteProblem::assemble(this->A, this->Dir, this->RHS, rhsonly);

  for (int i=0; i < this->get_num_dofs(); i++) RHS[i] += Dir[i];
}



