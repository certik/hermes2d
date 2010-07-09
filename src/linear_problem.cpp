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
LinearProblem::LinearProblem(WeakForm* wf_) : DiscreteProblem(wf_) {};
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

void LinearProblem::assemble(Matrix *A, Vector *RHS)
{
#ifdef H2D_COMPLEX
        this->assemble(A, RHS->get_c_array_cplx());
#else
        this->assemble(A, RHS->get_c_array());
#endif
}
void LinearProblem::assemble(Matrix *A, scalar *RHS)
{
    info("LinearProblem::assemble(Matrix *A, scalar *RHS)");
  // sanity checks
  int ndof = this->get_num_dofs();

  if (this->have_spaces == false)
    error("Before assemble(), you need to initialize spaces.");
  if (this->spaces == NULL) error("spaces = NULL in LinearProblem::assemble().");

  // Assemble the matrix A and vector RHS. If the problem is linear,
  // we need to subtract the vector Dir from RHS.
  scalar *Dir = new scalar[ndof];
  DiscreteProblem::assemble(A, Dir, RHS, false);

  for (int i=0; i < this->get_num_dofs(); i++) RHS[i] += Dir[i];
  delete Dir;
}

//// solve /////////////////////////////////////////////////////////////////////////////////////////

bool LinearProblem::solve(Matrix* mat, scalar* rhs, scalar* vec)
{
  int ndof = this->get_num_dofs();

  // sanity checks
  if (mat == NULL) error("matrix is NULL in LinearProblem::solve().");
  if (rhs == NULL) error("rhs is NULL in LinearProblem::solve().");
  if (vec == NULL) error("vec is NULL in LinearProblem::solve().");
  if (ndof == 0) error("ndof = 0 in LinearProblem::solve().");
  if (ndof != mat->get_size())
    error("Matrix size does not match ndof in in LinearProblem:solve().");

  // copy "rhs" into "vec", solve the matrix problem with "mat", "vec" ,
  // and save the result in "vec"
  memcpy(vec, rhs, sizeof(scalar) * ndof);
  bool flag = this->solve_matrix_problem(mat, vec);

  return flag;
}

bool LinearProblem::solve(Tuple<Solution*>sln_tuple)
{
  int n = sln_tuple.size();
  int ndof = this->get_num_dofs();

  // sanity checks
  if (n != this->wf->neq)
    error("Number of solutions does not match the number of equations in LinearProblem::solve().");
  if (this->Vec == NULL) error("Vec is NULL in LinearProblem::solve().");
  if (this->Vec_length != ndof || this->RHS_length != ndof || this->Dir_length != ndof)
    error("Length of vectors Vec, RHS or Dir does not match this->ndof in LinearProblem::solve().");
  if (ndof == 0) error("ndof = 0 in LinearProblem::solve().");
  if (ndof != this->A->get_size())
    error("Matrix size does not match vector length in LinearProblem:solve().");

  // solve the matrix problem with this->A and this->RHS and copy the 
  // result into this->Vec
  bool flag = this->solve(this->A, this->RHS, this->Vec);
  if (flag == false) return false; 

  // copy this->Vec into Solutions
  if (this->spaces == NULL) error("this->spaces == NULL in DiscreteProblem::solve().");
  if (this->pss == NULL) error("this->pss == NULL in DiscreteProblem::solve().");
  if (this->Vec == NULL) error("this->Vec == NULL in LinearProblem::solve().");
  for (int i = 0; i < n; i++)
  {
    if(this->spaces[i] == NULL) error("this->spaces[%d] == NULL in LinearProblem::solve().", i);
    if(this->spaces[i]->get_mesh() == NULL) error("this->spaces[%d]->get_mesh() == NULL in LinearProblem::solve().", i);
    sln_tuple[i]->set_fe_solution(this->spaces[i], this->pss[i], this->Vec);
    if(sln_tuple[i]->get_mesh() == NULL) error("sln_tuple[%d]->get_mesh() == NULL in LinearProblem::solve().\n", i);
  }

  return true;
}

