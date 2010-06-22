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

#ifndef __H2D_NONLINSYSTEM_H
#define __H2D_NONLINSYSTEM_H

#include "tuple.h"
#include "matrix_old.h"
#include "filter.h"
#include "views/scalar_view.h"
#include "views/order_view.h"
#include "function.h"

class Solution;
class MeshFunction;

H2D_API_USED_TEMPLATE(Tuple<Space*>); ///< Instantiated template. It is used to create a clean Windows DLL interface.
H2D_API_USED_TEMPLATE(Tuple<Solution*>); ///< Instantiated template. It is used to create a clean Windows DLL interface.

///
///
///
///
///
class H2D_API NonlinSystem : public LinSystem
{
public:

  /// Initializes the class and creates a zero initial coefficient vector.
  void init_nonlin();
  NonlinSystem();
  NonlinSystem(WeakForm* wf_, CommonSolver* solver_);
  NonlinSystem(WeakForm* wf_);                  // solver will be set to NULL and default solver will be used
  NonlinSystem(WeakForm* wf_, CommonSolver* solver_, Space* s_);
  NonlinSystem(WeakForm* wf_, Space* s_);        // solver will be set to NULL and default solver will be used
  NonlinSystem(WeakForm* wf_, CommonSolver* solver_, Tuple<Space*> spaces_);
  NonlinSystem(WeakForm* wf_, Tuple<Space*> spaces_);      // solver will be set to NULL and default solver will be used

  /// Frees the memory for the RHS, Dir and Vec vectors, and solver data.
  virtual void free();

  /// Adjusts the Newton iteration coefficient. The default value for alpha is 1.
  void set_alpha(double alpha) { this->alpha = alpha; }

  /// Assembles the Jacobian matrix and the residual vector (as opposed to 
  /// LibSystem::assemble() that constructs a stiffness matrix and load 
  /// vector). 
  void assemble(bool rhsonly = false);

  /// Performs one step of the Newton's method for an arbitrary number of equations, 
  /// stores the result in the given Solutions.
  bool solve(Tuple<Solution*> sln);
  bool solve(Solution* sln); // 1 equation
  bool solve(Solution* sln1, Solution* sln2); // 2 equations
  bool solve(Solution* sln1, Solution* sln2, Solution* sln3); // 3 equations

  /// Performs complete Newton's loop for a Tuple of solutions.
  bool solve_newton(Tuple<Solution*> u_prev, double newton_tol, int newton_max_iter,
                    bool verbose = false, Tuple<MeshFunction*> mesh_fns = Tuple<MeshFunction*>());

  /// Performs complete Newton's loop for one equation
  bool solve_newton(Solution* u_prev, double newton_tol, int newton_max_iter,
                    bool verbose = false, 
                    Tuple<MeshFunction*> mesh_fns = Tuple<MeshFunction*>())
  {
    return this->solve_newton(Tuple<Solution*>(u_prev), newton_tol, newton_max_iter, verbose, mesh_fns);
  }

  /// returns the L2-norm of the residual vector
  double get_residual_l2_norm() const { return res_l2; }

  /// returns the L1-norm of the residual vector
  double get_residual_l1_norm() const { return res_l1; }

  /// returns the L_inf-norm of the residual vector
  double get_residual_max_norm() const { return res_max; }
  
  virtual bool is_linear() { return false; }

protected:

  double alpha;
  double res_l2, res_l1, res_max;
};


#endif
