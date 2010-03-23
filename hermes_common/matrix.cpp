// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"

#include "python_api.h"

#define TINY 1e-20

/// Transposes an m by n matrix. If m != n, the array matrix in fact has to be
/// a square matrix of the size max(m, n) in order for the transpose to fit inside it.
template<typename T>
void transpose(T** matrix, int m, int n)
{
  int min = std::min(m, n);
  for (int i = 0; i < min; i++)
    for (int j = i+1; j < min; j++)
      std::swap(matrix[i][j], matrix[j][i]);

  if (m < n)
    for (int i = 0; i < m; i++)
      for (int j = m; j < n; j++)
        matrix[j][i] = matrix[i][j];
  else if (n < m)
    for (int i = n; i < m; i++)
      for (int j = 0; j < n; j++)
        matrix[j][i] = matrix[i][j];
}


/// Changes the sign of a matrix
template<typename T>
void chsgn(T** matrix, int m, int n)
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      matrix[i][j] = -matrix[i][j];
}


/// Given a matrix a[n][n], this routine replaces it by the LU decomposition of a rowwise
/// permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
/// indx[n] is an output vector that records the row permutation effected by the partial
/// pivoting; d is output as +-1 depending on whether the number of row interchanges was even
/// or odd, respectively. This routine is used in combination with lubksb to solve linear equations
/// or invert a matrix.
void ludcmp(double** a, int n, int* indx, double* d);

/// Solves the set of n linear equations AX = B. Here a[n][n] is input, not as the matrix
/// A but rather as its LU decomposition, determined by the routine ludcmp. indx[n] is input
/// as the permutation vector returned by ludcmp. b[n] is input as the right-hand side vector
/// B, and returns with the solution vector X. a, n, and indx are not modified by this routine
/// and can be left in place for successive calls with different right-hand sides b. This routine takes
/// into account the possibility that b will begin with many zero elements, so it is efficient for use
/// in matrix inversion.
//template<typename T>
void lubksb(double** a, int n, int* indx, double* b)
{
  int i, ip, j;
  double sum;

  for (i = 0; i < n; i++)
  {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    for (j = 0; j < i; j++) sum -= a[i][j]*b[j];
    b[i] = sum;
  }
  for (i = n-1; i >= 0; i--)
  {
    sum = b[i];
    for (j = i+1; j < n; j++) sum -= a[i][j]*b[j];
    b[i] = sum / a[i][i];
  }
}

/// Given a positive-definite symmetric matrix a[n][n], this routine constructs its Cholesky
/// decomposition, A = L*L^T . On input, only the upper triangle of a need be given; it is not
/// modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
/// elements which are returned in p[n].
void choldc(double **a, int n, double p[]);

/// Solves the set of n linear equations A*x = b, where a is a positive-definite symmetric matrix.
/// a[n][n] and p[n] are input as the output of the routine choldc. Only the lower
/// subdiagonal portion of a is accessed. b[n] is input as the right-hand side vector. The
/// solution vector is returned in x[n]. a, n, and p are not modified and can be left in place
/// for successive calls with different right-hand sides b. b is not modified unless you identify b and
/// x in the calling sequence, which is allowed. The right-hand side b can be complex, in which case
/// the solution x is also complex.
template<typename T>
void cholsl(double **a, int n, double p[], T b[], T x[])
{
  int i, k;
  T sum;

  for (i = 0; i < n; i++)
  {
    sum = b[i];
    k = i;
    while (--k >= 0)
      sum -= a[i][k] * x[k];
    x[i] = sum / p[i];
  }

  for (i = n-1; i >= 0; i--)
  {
    sum = x[i];
    k = i;
    while (++k < n)
      sum -= a[k][i] * x[k];
    x[i] = sum / p[i];
  }
}


void ludcmp(double** a, int n, int* indx, double* d)
{
  int i, imax = 0, j, k;
  double big, dum, sum, temp;
  double* vv = new double[n];

  *d = 1.0;
  for (i = 0; i < n; i++)
  {
    big=0.0;
    for (j = 0; j < n; j++) 
      if ((temp = fabs(a[i][j])) > big) 
        big = temp;
    if (big == 0.0) _error("Singular matrix!");
    vv[i] = 1.0 / big;
  }
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < j; i++)
    {
      sum = a[i][j];
      for (k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++)
    {
      sum = a[i][j];
      for (k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if ((dum = vv[i]*fabs(sum)) >= big)
      {
        big = dum;
        imax = i;
      }
    }
    if (j != imax)
    {
      for (k = 0; k < n; k++)
      {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0) a[j][j] = TINY;
    if (j != n-1) 
    {
      dum = 1.0 / (a[j][j]);
      for (i = j+1; i < n; i++) a[i][j] *= dum;
    }
  }
  delete [] vv;
}


void choldc(double **a, int n, double p[])
{
  int i, j, k;
  for (i = 0; i < n; i++)
  {
    for (j = i; j < n; j++)
    {
      double sum = a[i][j];
      k = i;
      while (--k >= 0)
        sum -= a[i][k] * a[j][k];
      
      if (i == j)
      {
        if (sum <= 0.0)
          _error("CHOLDC failed!");
        else
          p[i] = sqrt(sum);
      }
      else
        a[j][i] = sum / p[i];
    }
  }
}

Matrix::Matrix()
{
    this->p = new Python();
}

Matrix::~Matrix()
{
    delete this->p;
}

CSRMatrix::CSRMatrix(Matrix *m):Matrix() {
    if (dynamic_cast<CooMatrix*>(m))
        this->add_from_CooMatrix((CooMatrix*)m);
    else if (dynamic_cast<CSCMatrix*>(m))
        this->add_from_CSCMatrix((CSCMatrix*)m);
    else if (dynamic_cast<DenseMatrix*>(m))
        this->add_from_DenseMatrix((DenseMatrix*)m);
    else
        _error("Matrix type not supported.");
}

void CSRMatrix::add_from_CooMatrix(CooMatrix *m)
{
    p->push("m", c2py_CooMatrix(m));
    p->exec("n = m.to_scipy_coo().tocsr()");
    p->exec("A = n.data");
    p->exec("IA = n.indptr");
    p->exec("JA = n.indices");
    //Note: these arrays are allocated by numpy and the pointers A, IA, JA
    //point to them (inplace) so we will leave the deallocation to numpy
    this->_is_complex = m->is_complex();
    if (this->is_complex())
        numpy2c_double_complex_inplace(p->pull("A"), &(this->A_cplx),
                &(this->nnz));
    else
        numpy2c_double_inplace(p->pull("A"), &(this->A), &(this->nnz));
    numpy2c_int_inplace(p->pull("IA"), &(this->IA), &(this->size));
    numpy2c_int_inplace(p->pull("JA"), &(this->JA), &(this->nnz));
    this->size--;
    this->deallocate_arrays = false;
}

void CSRMatrix::add_from_CSCMatrix(CSCMatrix *m)
{
    p->push("m", c2py_CSCMatrix(m));
    p->exec("n = m.to_scipy_csc().tocsr()");
    p->exec("A = n.data");
    p->exec("IA = n.indptr");
    p->exec("JA = n.indices");
    //Note: these arrays are allocated by numpy and the pointers A, IA, JA
    //point to them (inplace) so we will leave the deallocation to numpy
    this->_is_complex = m->is_complex();
    if (this->is_complex())
        numpy2c_double_complex_inplace(p->pull("A"), &(this->A_cplx),
                &(this->nnz));
    else
        numpy2c_double_inplace(p->pull("A"), &(this->A), &(this->nnz));
    numpy2c_int_inplace(p->pull("IA"), &(this->IA), &(this->size));
    numpy2c_int_inplace(p->pull("JA"), &(this->JA), &(this->nnz));
    this->size--;
    this->deallocate_arrays = false;
}

void CSRMatrix::print()
{
    p->push("m", c2py_CSRMatrix(this));
    p->exec("S = str(m.to_scipy_csr())");
    printf("%s\n", py2c_str(p->pull("S")));
}

void CSCMatrix::add_from_CooMatrix(CooMatrix *m)
{
    p->push("m", c2py_CooMatrix(m));
    p->exec("n = m.to_scipy_coo().tocsc()");
    p->exec("A = n.data");
    p->exec("IA = n.indices");
    p->exec("JA = n.indptr");
    //Note: these arrays are allocated by numpy and the pointers A, IA, JA
    //point to them (inplace) so we will leave the deallocation to numpy
    this->_is_complex = m->is_complex();
    if (this->is_complex())
        numpy2c_double_complex_inplace(p->pull("A"), &(this->A_cplx),
                &(this->nnz));
    else
        numpy2c_double_inplace(p->pull("A"), &(this->A), &(this->nnz));
    numpy2c_int_inplace(p->pull("IA"), &(this->IA), &(this->nnz));
    numpy2c_int_inplace(p->pull("JA"), &(this->JA), &(this->size));
    this->size--;
}

void CSCMatrix::add_from_CSRMatrix(CSRMatrix *m)
{
    p->push("m", c2py_CSRMatrix(m));
    p->exec("n = m.to_scipy_csr().tocsc()");
    p->exec("A = n.data");
    p->exec("IA = n.indices");
    p->exec("JA = n.indptr");
    //Note: these arrays are allocated by numpy and the pointers A, IA, JA
    //point to them (inplace) so we will leave the deallocation to numpy
    this->_is_complex = m->is_complex();
    if (this->is_complex())
        numpy2c_double_complex_inplace(p->pull("A"), &(this->A_cplx),
                &(this->nnz));
    else
        numpy2c_double_inplace(p->pull("A"), &(this->A), &(this->nnz));
    numpy2c_int_inplace(p->pull("IA"), &(this->IA), &(this->nnz));
    numpy2c_int_inplace(p->pull("JA"), &(this->JA), &(this->size));
    this->size--;
}

void CSCMatrix::print()
{
    p->push("m", c2py_CSCMatrix(this));
    p->exec("S = str(m.to_scipy_csc())");
    printf("%s\n", py2c_str(p->pull("S")));
}

void CooMatrix::print()
{
    p->push("m", c2py_CooMatrix(this));
    p->exec("S = str(m.to_scipy_coo())");
    printf("%s\n", py2c_str(p->pull("S")));
}

template <> Triple<cplx> *CooMatrix::get_list<cplx>() {
    return this->list_cplx;
}

void CooMatrix::set_zero() {
    if (this->_is_complex) {
        this->free_data<cplx>();
        this->list_cplx = NULL;
        this->list_last_cplx = NULL;
    } else {
        this->free_data<double>();
        this->list = NULL;
        this->list_last = NULL;
    }
    this->size = 0;
}


// matrix vector multiplication
void mat_dot(Matrix* A, double* x, double* result, int n_dof)
{
  A->times_vector(x, result, n_dof);
}

// vector vector multiplication
double vec_dot(double* r, double* s, int n_dof)
{
  double result = 0;
  for (int i=0; i < n_dof; i++) result += r[i]*s[i];
  return result;
}


void solve_linear_system_dense_lu(Matrix *mat, double *res)
{
    DenseMatrix *dmat = dynamic_cast<DenseMatrix*>(mat);
    bool free_dmat = false;
    if (dmat == NULL) {
        dmat = new DenseMatrix(mat);
        free_dmat = true;
    }
    int n = dmat->get_size();
    int *indx = new int[n];
    double **_mat = dmat->get_mat();
    double d;
    ludcmp(_mat, n, indx, &d);
    lubksb(_mat, n, indx, res);
    if (free_dmat)
        delete dmat;
}

// Standard CG method starting from zero vector
// (because we solve for the increment)
// x... comes as right-hand side, leaves as solution
int solve_linear_system_cg(Matrix* A, double *x,
                           double matrix_solver_tol,
                           int matrix_solver_maxiter)
{
  int n_dof = A->get_size();
  double *r = new double[n_dof];
  double *p = new double[n_dof];
  double *help_vec = new double[n_dof];
  if (r == NULL || p == NULL || help_vec == NULL) {
    _error("a vector could not be allocated in solve_linear_system_iter().");
  }
  // r = b - A*x0  (where b is x and x0 = 0)
  for (int i=0; i < n_dof; i++) r[i] = x[i];
  // p = r
  for (int i=0; i < n_dof; i++) p[i] = r[i];

  // setting initial condition x = 0
  for (int i=0; i < n_dof; i++) x[i] = 0;

  // CG iteration
  int iter_current = 0;
  double tol_current;
  while (1) {
    mat_dot(A, p, help_vec, n_dof);
    double r_times_r = vec_dot(r, r, n_dof);
    double alpha = r_times_r / vec_dot(p, help_vec, n_dof); 
    for (int i=0; i < n_dof; i++) {
      x[i] += alpha*p[i];
      r[i] -= alpha*help_vec[i];
    }
    double r_times_r_new = vec_dot(r, r, n_dof);
    iter_current++;
    tol_current = sqrt(r_times_r_new);
    if (tol_current < matrix_solver_tol 
        || iter_current >= matrix_solver_maxiter) break;
    double beta = r_times_r_new/r_times_r;
    for (int i=0; i < n_dof; i++) p[i] = r[i] + beta*p[i];
  }
  int flag;
  if (tol_current <= matrix_solver_tol) flag = 1;
  else flag = 0;
  if (r != NULL) delete [] r;
  if (p != NULL) delete [] p;
  if (help_vec != NULL) delete [] help_vec;

  printf("CG (regular) made %d iteration(s) (tol = %g)\n", 
         iter_current, tol_current);

  return flag;
}
