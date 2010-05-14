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

#ifndef __H2D_SOLVER_UMFPACK_H
#define __H2D_SOLVER_UMFPACK_H

//// UmfpackSolver /////////////////////////////////////////////////////////////////////////////////

#include "common.h"
#include "solver.h"

/* This is needed for backwards compatibility only. This class does nothing and
 * should not be used. Use the solvers via hermes_common.
 */

class UmfpackSolver : public Solver {
public:
    UmfpackSolver() {}

protected:
    virtual bool is_row_oriented() {}
    virtual bool handles_symmetry() {}
    virtual bool solve(void*, int, int*, int*, scalar*, bool, scalar*,
            scalar*){}
};

#endif
