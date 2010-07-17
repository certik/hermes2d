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

void LinearProblem::assemble(Matrix* mat_ext, Vector* rhs_ext, bool rhsonly, bool is_complex)
{
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof == 0 in LinearProblem::assemble().");
  Vector* dir_ext = new AVector(ndof, is_complex);
  // the vector dir represents the contribution of the Dirichlet lift, 
  // and for linear problems  it has to be subtracted from the right hand side
  DiscreteProblem::assemble(mat_ext, dir_ext, rhs_ext, rhsonly);
  // FIXME: Do we really need to handle the real and complex cases separately?
  if (is_complex) for (int i=0; i < ndof; i++) rhs_ext->add(i, -dir_ext->get_cplx(i));
  else for (int i=0; i < ndof; i++) rhs_ext->add(i, -dir_ext->get(i));
  delete dir_ext;
}

// FIXME: We need to unify the type for Python and 
// C++ solvers. Right now Solver and CommonSolver
// are incompatible.
void init_matrix_solver(MatrixSolverType matrix_solver, int ndof, 
                        Matrix* &mat, Vector* &rhs, CommonSolver* &solver, bool is_complex) 
{
  // Initialize stiffness matrix, load vector, and matrix solver.
  // UMFpack.
  CooMatrix* mat_umfpack = new CooMatrix(ndof, is_complex);
  Vector* rhs_umfpack = new AVector(ndof, is_complex);
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

// Solve a typical linear problem (without automatic adaptivity).
// Feel free to adjust this function for more advanced applications.
bool solve_linear(Tuple<Space *> spaces, WeakForm* wf, Tuple<Solution *> solutions, 
                  MatrixSolverType matrix_solver, bool is_complex) 
{
  // Initialize the linear problem.
  LinearProblem lp(wf, spaces);
  //info("ndof = %d", lp.get_num_dofs());

  // Select matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;
  init_matrix_solver(matrix_solver, lp.get_num_dofs(), mat, rhs, solver, is_complex);

  // Assemble stiffness matrix and rhs.
  bool rhsonly = false;
  lp.assemble(mat, rhs, rhsonly, is_complex);

  //mat->print();

  // Solve the matrix problem.
  if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

  // Convert coefficient vector into a Solution.
  for (int i=0; i<solutions.size(); i++) {
    solutions[i]->set_fe_solution(spaces[i], rhs);
  }
}

int get_num_dofs(Tuple<Space *> spaces) 
{
  int ndof = 0;
  for (int i=0; i<spaces.size(); i++) {
    ndof += spaces[i]->get_num_dofs();
  }
  return ndof;
}

// Solve a typical linear problem (without automatic adaptivity).
// Feel free to adjust this function for more advanced applications.
bool solve_linear_adapt(Tuple<Space *> spaces, WeakForm* wf, Tuple<Solution *> solutions, 
                        MatrixSolverType matrix_solver, Tuple<int> proj_norms, 
                        RefinementSelectors::Selector* selector, AdaptivityParamType* apt,
                        const int sln_win_geom[4], const int mesh_win_geom[4], 
                        bool is_complex) 
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Adaptivity parameters.
  double err_stop = apt->err_stop; 
  int ndof_stop = apt->ndof_stop;
  double threshold = apt->threshold;
  int strategy = apt->strategy; 
  int mesh_regularity = apt->mesh_regularity;
  double to_be_processed = apt->to_be_processed;

  // Number of physical fields in the problem.
  int num_comps = spaces.size();

  // Calculate the number of degreeso of freedom 
  int ndof = get_num_dofs(spaces);

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  init_matrix_solver(matrix_solver, ndof, mat, rhs, solver, is_complex);

  // Initialize views.
  ScalarView sview("Solution", sln_win_geom);
  OrderView  oview("Mesh", mesh_win_geom);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu;

  // Declare coarse mesh and reference solutions
  Tuple<Solution *> ref_slns;
  for (int i = 0; i < num_comps; i++) {
    Solution *ref_sln = new Solution();
    ref_slns.push_back(ref_sln);
  } 

  // Conversion from Tuple<Solution *> to Tuple<MeshFunction *>
  // so that project_global() below compiles. 
  Tuple<MeshFunction *> ref_slns_mf;
  for (int i = 0; i < num_comps; i++) {
    MeshFunction *s = (MeshFunction*)ref_slns[i];
    ref_slns_mf.push_back(s);
  }

  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
    if (num_comps == 1) info("Solving on reference mesh.");
    else info("Solving on reference meshes.");

    // Construct globally refined reference mesh(es)
    // and setup reference space(s).
    Tuple<Space *> ref_spaces;
    for (int i = 0; i < num_comps; i++) {
      Mesh *ref_mesh = new Mesh();
      ref_mesh->copy(spaces[i]->get_mesh());
      ref_mesh->refine_all_elements();
      ref_spaces.push_back(spaces[i]->dup(ref_mesh));
      int order_increase = 1;
      ref_spaces[i]->copy_orders(spaces[i], order_increase);
    }

    // Solve the reference problem.
    solve_linear(ref_spaces, wf, ref_slns, matrix_solver);

    // Project the reference solution on the coarse mesh.
    if (num_comps == 1) info("Projecting reference solution on coarse mesh.");
    else info("Projecting reference solutions on coarse meshes.");
    project_global(spaces, ref_slns_mf, solutions, proj_norms, is_complex);

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution (first component only).
    sview.show(solutions[0]);
    oview.show(spaces[0]);

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors and total error estimate.
    info("Calculating error.");
    H1Adapt hp(spaces);
    hp.set_solutions(solutions, ref_slns);
    double err_est = hp.calc_error() * 100;

    // Report results.
    info("ndof: %d, ref_ndof: %d, err_est: %g%%", 
         get_num_dofs(spaces), get_num_dofs(ref_spaces), err_est);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(get_num_dofs(spaces), err_est);
    graph_dof.save("conv_dof.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu.dat");

    // If err_est too large, adapt the mesh.
    if (err_est < err_stop) done = true;
    else {
      if (num_comps == 1) info("Adapting the coarse mesh.");
      else info("Adapting the coarse meshes.");
      done = hp.adapt(selector, threshold, strategy, mesh_regularity, to_be_processed);

      if (get_num_dofs(spaces) >= ndof_stop) done = true;
    }

    as++;
  }
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());
}
