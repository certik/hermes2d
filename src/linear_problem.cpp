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
#include "tuple.h"

LinearProblem::LinearProblem() : DiscreteProblem() {};
LinearProblem::LinearProblem(WeakForm* wf_) : DiscreteProblem(wf_) {};
LinearProblem::LinearProblem(WeakForm* wf_, Space* s_) : DiscreteProblem(wf_, s_) {};
LinearProblem::LinearProblem(WeakForm* wf_, Tuple<Space*> spaces_) : DiscreteProblem(wf_, spaces_) {};
LinearProblem::~LinearProblem() {};

void LinearProblem::assemble(Matrix* mat_ext, Vector* rhs_ext, bool rhsonly)
{
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof == 0 in LinearProblem::assemble().");
  Vector* dir_ext = new AVector(ndof);
  // the vector dir represents the contribution of the Dirichlet lift, 
  // and it has to be subtracted from the right hand side for linear problems
  DiscreteProblem::assemble(mat_ext, dir_ext, rhs_ext, rhsonly);
  for (int i=0; i < ndof; i++) rhs_ext->add(i, -dir_ext->get(i));
  delete dir_ext;
}

// FIXME: We need to unify the type for Python and 
// C++ solvers. Right now Solver and CommonSolver
// are incompatible.
void init_matrix_solver(MatrixSolverType matrix_solver, int ndof, 
                        Matrix* &mat, Vector* &rhs, CommonSolver* &solver) 
{
  // Initialize stiffness matrix, load vector, and matrix solver.
  // UMFpack.
  CooMatrix* mat_umfpack = new CooMatrix(ndof);
  Vector* rhs_umfpack = new AVector(ndof);
  CommonSolverSciPyUmfpack* solver_umfpack = new CommonSolverSciPyUmfpack();
    //CommonSolverSciPyUmfpack* solver_umfpack = new CommonSolverSciPyUmfpack();
  // PETSc.
  /* FIXME - PETSc solver needs to be ported from H3D.
  PetscMatrix mat_petsc(ndof);
  PetscVector rhs_petsc(ndof);
  PetscLinearSolver solver_petsc;
  */
  // MUMPS. 
  // FIXME - PETSc solver needs to be ported from H3D.
  /*
  MumpsMatrix mat_mumps(ndof);
  MumpsVector rhs_mumps(ndof);
  MumpsSolver solver_mumps;
  */
  
  switch (matrix_solver) {
    case SOLVER_UMFPACK: 
      mat = mat_umfpack;
      rhs = rhs_umfpack;
      solver = solver_umfpack;
      break;
    case SOLVER_PETSC:  
      error("Petsc solver not implemented yet.");
      /*
      mat = &mat_petsc;
      rhs = &rhs_petsc;
      solver = &solver_petsc;
      */
      break;
    case SOLVER_MUMPS:  
      error("MUMPS solver not implemented yet.");
      /*
      mat = &mat_mumps;
      rhs = &rhs_mumps;
      solver = &solver_mumps;
      */
      break;
    default: error("Bad matrix solver in init_matrix_solver().");
  }
}

// Shortcut to solve linear problems.
bool solve_linear(Tuple<Space *> spaces, WeakForm* wf, Tuple<Solution *> solutions, 
                  MatrixSolverType matrix_solver) 
{
  // Initialize the linear problem.
  LinearProblem lp(wf, spaces);
  info("ndof = %d", lp.get_num_dofs());

  // Select matrix solver.
  Matrix* mat;
  Vector* rhs;
  CommonSolver* solver;  // FIXME: this should be just Solver, same for
                         // Python and C++ solvers. 
  init_matrix_solver(matrix_solver, lp.get_num_dofs(), mat, rhs, solver);

  // Assemble stiffness matrix and rhs.
  lp.assemble(mat, rhs);

  //mat->print();

  // Solve the matrix problem.
  if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

  // Convert coefficient vector into a Solution.
  for (int i=0; i<solutions.size(); i++) {
    solutions[i]->set_fe_solution(spaces[i], rhs);
  }
}
