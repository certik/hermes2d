// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES_COMMON_SOLVERS_H
#define __HERMES_COMMON_SOLVERS_H

class Matrix;
class Vector;

// abstract class
class CommonSolver
{
public:
    virtual bool solve(Matrix* mat, Vector* res) = 0;
    inline char *get_log() { return log; }

private:
    char *log;
};

// c++ cg
class CommonSolverCG : public CommonSolver
{
public:
    bool solve(Matrix* mat, Vector* res, double tol = 1e-6, int maxiter = 1000);
    bool solve(Matrix* mat, Vector* res)
    {
      this->solve(mat, res, 1e-6, 1000);
    };
};
inline bool solve_linear_system_cg(Matrix *mat, Vector *res,
                                   double tolerance = 1e-6,
                                   int maxiter = 1000)
{
    CommonSolverCG solver;
    return solver.solve(mat, res, tolerance, maxiter);
}

// c++ lu
class CommonSolverDenseLU : public CommonSolver
{
public:
    bool solve(Matrix *mat, Vector *res);
};
inline void solve_linear_system_dense_lu(Matrix *mat, Vector *res)
{
    CommonSolverDenseLU solver;
    solver.solve(mat, res);
}

// c++ umfpack - optional
class CommonSolverUmfpack : public CommonSolver
{
public:
    bool solve(Matrix *mat, Vector *res);
    bool solve_real(Matrix *mat, double *res);
    bool solve_cplx(Matrix *mat, cplx *res);
};
inline void solve_linear_system_umfpack(Matrix *mat, Vector *res)
{
    CommonSolverUmfpack solver;
    solver.solve(mat, res);
}

// c++ sparselib - optional
class CommonSolverSparseLib : public CommonSolver
{
public:
    enum CommonSolverSparseLibSolver
    {
        CommonSolverSparseLibSolver_ConjugateGradientSquared,
        CommonSolverSparseLibSolver_RichardsonIterativeRefinement
    };

    CommonSolverSparseLib()
    {
        tolerance = 1e-8;
        maxiter = 1000;
        method = CommonSolverSparseLibSolver_ConjugateGradientSquared;
    }

    bool solve(Matrix *mat, Vector *res);
    inline void set_tolerance(double tolerance) { this->tolerance = tolerance; }
    inline void set_maxiter(int maxiter) { this->maxiter = maxiter; }
    inline void set_method(CommonSolverSparseLibSolver method) { this->method = method; }

private:
    double tolerance;
    int maxiter;
    CommonSolverSparseLibSolver method;
};
inline void solve_linear_system_sparselib_cgs(Matrix *mat, Vector *res, double tolerance = 1e-8, int maxiter = 1000)
{
    CommonSolverSparseLib solver;
    solver.set_method(CommonSolverSparseLib::CommonSolverSparseLibSolver_ConjugateGradientSquared);
    solver.set_tolerance(tolerance);
    solver.set_maxiter(maxiter);
    solver.solve(mat, res);
}
inline void solve_linear_system_sparselib_ir(Matrix *mat, Vector *res, double tolerance = 1e-8, int maxiter = 1000)
{
    CommonSolverSparseLib solver;
    solver.set_method(CommonSolverSparseLib::CommonSolverSparseLibSolver_RichardsonIterativeRefinement);
    solver.set_tolerance(tolerance);
    solver.set_maxiter(maxiter);
    solver.solve(mat, res);
}

// c++ superlu - optional
class CommonSolverSuperLU : public CommonSolver
{
public:
    bool solve(Matrix *mat, Vector *res);
    bool solve2(Matrix *mat, Vector *res);
};
inline void solve_linear_system_superlu(Matrix *mat, Vector *res)
{
    CommonSolverSuperLU solver;
    solver.solve(mat, res);
}

// python numpy - optional
class CommonSolverNumPy : public CommonSolver
{
public:
    bool solve(Matrix *mat, Vector *res);
};
inline void solve_linear_system_numpy(Matrix *mat, Vector *res)
{
    CommonSolverNumPy solver;
    solver.solve(mat, res);
}

// python scipy - optional
class CommonSolverSciPyUmfpack : public CommonSolver
{
public:
    bool solve(Matrix *mat, Vector *res);
};
inline void solve_linear_system_scipy_umfpack(Matrix *mat, Vector *res)
{
    CommonSolverSciPyUmfpack solver;
    solver.solve(mat, res);
}

class CommonSolverSciPyCG : public CommonSolver
{
public:
    bool solve(Matrix *mat, Vector *res);
};
inline void solve_linear_system_scipy_cg(Matrix *mat, Vector *res)
{
    CommonSolverSciPyCG solver;
    solver.solve(mat, res);
}

class CommonSolverSciPyGMRES : public CommonSolver
{
public:
    bool solve(Matrix *mat, Vector *res);
};
inline void solve_linear_system_scipy_gmres(Matrix *mat, Vector *res)
{
    CommonSolverSciPyGMRES solver;
    solver.solve(mat, res);
}

#endif
