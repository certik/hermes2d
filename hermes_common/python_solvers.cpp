// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"
#include "solvers.h"

#ifdef COMMON_WITH_SCIPY
#include "python_api.h"

bool CommonSolverNumPy::solve(Matrix *mat, Vector *res)
{
  //printf("NumPy solver\n");

    double* vec_double;
    cplx* vec_cplx;
    if (sizeof(scalar) == sizeof(double)) vec_double = (double*)res->get_c_array();
    else vec_cplx = (cplx*)res->get_c_array();

    CSRMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSRMatrix(&M));
    if (sizeof(scalar) == sizeof(double)) p->push("rhs", c2numpy_double_inplace(vec_double, mat->get_size()));
    else p->push("rhs", c2numpy_double_complex_inplace(vec_cplx, mat->get_size()));
    p->exec("A = m.to_scipy_csr().todense()");
    p->exec("from numpy.linalg import solve");
    p->exec("x = solve(A, rhs)");
    int n;
    if (sizeof(scalar) == sizeof(double)) {
      double *x_double;
      numpy2c_double_inplace(p->pull("x"), &x_double, &n);
      memcpy(vec_double, x_double, n*sizeof(scalar));
    }
    else {
      cplx *x_cplx; 
      numpy2c_double_complex_inplace(p->pull("x"), &x_cplx, &n);
      memcpy(vec_cplx, x_cplx, n*sizeof(scalar));
    }
    delete p;
}

bool CommonSolverSciPyUmfpack::solve(Matrix *mat, Vector *res)
{
  //printf("SciPy UMFPACK solver\n");

    double* vec_double;
    cplx* vec_cplx;
    if (sizeof(scalar) == sizeof(double)) vec_double = (double*)res->get_c_array();
    else vec_cplx = (cplx*)res->get_c_array();

    CSCMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSCMatrix(&M));
    if (sizeof(scalar) == sizeof(double)) p->push("rhs", c2numpy_double_inplace(vec_double, mat->get_size()));
    else p->push("rhs", c2numpy_double_complex_inplace(vec_cplx, mat->get_size()));
    p->exec("A = m.to_scipy_csc()");
    p->exec("from scipy.sparse.linalg import spsolve");
    // Turn off warnings in spsolve (only there)
    p->exec("import warnings");
    p->exec("with warnings.catch_warnings():\n"
            "    warnings.simplefilter('ignore')\n"
            "    x = spsolve(A, rhs)");   
    int n;
    if (sizeof(scalar) == sizeof(double)) {
      double *x_double;
      numpy2c_double_inplace(p->pull("x"), &x_double, &n);
      memcpy(vec_double, x_double, n*sizeof(scalar));
    }
    else {
      cplx *x_cplx;
      numpy2c_double_complex_inplace(p->pull("x"), &x_cplx, &n);
      memcpy(vec_cplx, x_cplx, n*sizeof(scalar));

    }
    delete p;
    return true;
}

bool CommonSolverSciPyCG::solve(Matrix *mat, Vector *res)
{
  //printf("SciPy CG solver\n");

    double* vec_double;
    cplx* vec_cplx;
    if (sizeof(scalar) == sizeof(double)) vec_double = (double*)res->get_c_array();
    else vec_cplx = (cplx*)res->get_c_array();

    CSRMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSRMatrix(&M));
    if (sizeof(scalar) == sizeof(double)) p->push("rhs", c2numpy_double_inplace(vec_double, mat->get_size()));
    else p->push("rhs", c2numpy_double_complex_inplace(vec_cplx, mat->get_size()));
    p->exec("A = m.to_scipy_csr()");
    p->exec("from scipy.sparse.linalg import cg");
    p->exec("x, res = cg(A, rhs)");
    int n;
    if (sizeof(scalar) == sizeof(double)) {
      double *x_double;
      numpy2c_double_inplace(p->pull("x"), &x_double, &n);
      memcpy(vec_double, x_double, n*sizeof(scalar));
    }
    else {
      cplx *x_cplx; 
      numpy2c_double_complex_inplace(p->pull("x"), &x_cplx, &n);
      memcpy(vec_cplx, x_cplx, n*sizeof(scalar));
    }
    delete p;
}

bool CommonSolverSciPyGMRES::solve(Matrix *mat, Vector *res)
{
  //printf("SciPy GMRES solver\n");

    double* vec_double;
    cplx* vec_cplx;
    if (sizeof(scalar) == sizeof(double)) vec_double = (double*)res->get_c_array();
    else vec_cplx = (cplx*)res->get_c_array();

    CSRMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSRMatrix(&M));
    if (sizeof(scalar) == sizeof(double)) p->push("rhs", c2numpy_double_inplace(vec_double, mat->get_size()));
    else p->push("rhs", c2numpy_double_complex_inplace(vec_cplx, mat->get_size()));
    p->exec("A = m.to_scipy_csr()");
    p->exec("from scipy.sparse.linalg import gmres");
    p->exec("x, res = gmres(A, rhs)");
    int n;
    if (sizeof(scalar) == sizeof(double)) {
      double *x_double;
      numpy2c_double_inplace(p->pull("x"), &x_double, &n);
      memcpy(vec_double, x_double, n*sizeof(scalar));
    }
    else {
      cplx *x_cplx; 
      numpy2c_double_complex_inplace(p->pull("x"), &x_cplx, &n);
      memcpy(vec_cplx, x_cplx, n*sizeof(scalar));
    }
    delete p;
}

#else

bool CommonSolverNumPy::solve(Matrix *mat, double *res)
{
    _error("CommonSolverNumPy::solve(Matrix *mat, double *res) not implemented.");
}

bool CommonSolverNumPy::solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverNumPy::solve(Matrix *mat, cplx *res) not implemented.");
}

bool CommonSolverSciPyUmfpack::solve(Matrix *mat, double *res)
{
    _error("CommonSolverSciPyUmfpack::solve(Matrix *mat, double *res) not implemented.");
}

bool CommonSolverSciPyUmfpack::solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSciPyUmfpack::solve(Matrix *mat, cplx *res) not implemented.");
}

bool CommonSolverSciPyCG::solve(Matrix *mat, double *res)
{
    _error("CommonSolverSciPyCG::solve(Matrix *mat, double *res) not implemented.");
}

bool CommonSolverSciPyCG::solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSciPyCG::solve(Matrix *mat, cplx *res) not implemented.");
}

bool CommonSolverSciPyGMRES::solve(Matrix *mat, double *res)
{
    _error("CommonSolverSciPyGMRES::solve(Matrix *mat, double *res) not implemented.");
}

bool CommonSolverSciPyGMRES::solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSciPyGMRES::solve(Matrix *mat, cplx *res) not implemented.");
}


#endif
