// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES_COMMON_MATRIX_H
#define __HERMES_COMMON_MATRIX_H

#include <math.h>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <math.h>
#include <string.h>
#include <complex>
#include <map>

#include "scalar.h"

class Matrix;

#include "solvers.h"

// printf debug information about the stiffness/Jacobi matrix
#define _error(x) throw std::runtime_error(x)

class CooMatrix;
class CSRMatrix;
class CSCMatrix;

/// Creates a new (full) matrix with m rows and n columns with entries of the type T.
/// The entries can be accessed by matrix[i][j]. To delete the matrix, just
/// do "delete matrix".
template<typename T>
T** _new_matrix(int m, int n = 0)
{
    if (!n) n = m;
    T** vec = (T**) new char[sizeof(T*)*m + sizeof(T)*m*n];
    if (vec == NULL) _error("Out of memory.");
    T* row = (T*) (vec + m);
    for (int i = 0; i < m; i++, row += n)
        vec[i] = row;
    return vec;
}

class Matrix {
public:
    Matrix() {}
    virtual ~Matrix() {}

    inline virtual void init() { free_data(); }
    virtual void free_data() = 0;

    virtual void set_zero() = 0;

    inline virtual int get_size() { return this->size; }
    virtual void print() = 0;

    virtual void add(int m, int n, scalar v) = 0;

    virtual void add_block(int *iidx, int ilen, int *jidx, int jlen, scalar** mat)
    {
        for (int i = 0; i < ilen; i++)
            for (int j=0; j < jlen; j++)
                if (iidx[i] >= 0 && jidx[j] >= 0)
                    this->add(iidx[i], jidx[j], mat[i][j]);
    }
    virtual scalar get(int m, int n) = 0;
    virtual void copy_into(Matrix *m) = 0;

    virtual void times_vector(scalar* vec, scalar* result, int rank)
    {
        _error("internal error: times_vector() not implemented.");
    }

protected:
    int size;
};

class Vector {
public:
    Vector() {}
    virtual ~Vector() {}
    inline virtual int get_size() { return this->size; }
    virtual void print() = 0;

    virtual void add(int m, scalar v) = 0;
    virtual void add_block(int *iidx, int ilen, scalar* vec)
    {
        for (int i = 0; i < ilen; i++)
            if (iidx[i] >= 0)
                this->add(iidx[i], vec[i]);
    }
    virtual scalar set(int m, scalar val) = 0;
    virtual scalar get(int m) = 0;
    virtual scalar* get_c_array()
    {
        _error("internal error: get_c_array() not implemented.");
    }
    virtual void realloc_and_erase(int new_length)
    {
        _error("internal error: realloc_and_erase() not implemented.");
    };
protected:
    int size;
};

// print vector - int
void print_vector(const char *label, int *value, int size);
// print vector - double
void print_vector(const char *label, double *value, int size);
// print vector - complex
void print_vector(const char *label, cplx *value, int size);


// Uses a C++ array as the internal implementation
class AVector: public Vector {
public:
    AVector(int n) {
        this->size = n;
        this->v = new scalar[n];
        for (int i=0; i < n; i++) this->v[i] = 0;
    }
    virtual ~AVector() {
        delete[] this->v;
    }
    virtual void print() {
        print_vector("", this->v, this->get_size());
    }

    virtual void add(int m, scalar v) {
        this->v[m] += v;
    }
    virtual scalar set(int m, scalar val) {
        this->v[m] = val;
    }
    virtual scalar get(int m) {
        return this->v[m];
    }
    virtual scalar *get_c_array()
    {
        return this->v;
    }
    virtual void realloc_and_erase(int new_length)
    {
      this->v = (scalar*)realloc(this->v, new_length*sizeof(scalar));
      memset(this->v, 0, new_length*sizeof(scalar));
      this->size = new_length;
    }
;
private:
    scalar *v;
};

// **********************************************************************************************************

class CooMatrix : public Matrix {
public:
    CooMatrix();
    CooMatrix(int size);
    CooMatrix(Matrix *m);
    CooMatrix(CooMatrix *m);
    CooMatrix(CSRMatrix *m);
    CooMatrix(CSCMatrix *m);
    ~CooMatrix();

    inline virtual void init() { free_data(); }
    virtual void free_data();

    virtual void set_zero()
    {
        _error("CooMatrix::set_zero() not implemented.");
    }

    virtual int get_nnz();
    virtual void print();

    void add_from_csr(CSRMatrix *m);
    void add_from_csc(CSCMatrix *m);

    virtual void add(int m, int n, scalar v);
    void get_row_col_data(int *row, int *col, scalar *data);
    void get_row_col_data(int *row, int *col, double *data_real, double *data_imag);

    virtual void copy_into(Matrix *m);

    inline virtual scalar get(int m, int n) { return A[m][n]; }

    virtual void times_vector(scalar* vec, scalar* result, int rank);

protected:
    std::map<size_t, std::map<size_t, scalar> > A;
};

// **********************************************************************************************************

class DenseMatrix : public Matrix
{
public:
    DenseMatrix(Matrix *m);
    DenseMatrix(CooMatrix *m);
    DenseMatrix(int size);
    ~DenseMatrix();

    virtual void init();

    virtual void free_data();
    virtual void set_zero();

    inline virtual void add(int m, int n, scalar v) { this->A[m][n] += v; }

    void add_from_coo(CooMatrix *m);

    inline virtual scalar get(int m, int n) { return this->A[m][n]; }

    virtual void copy_into(Matrix *m)
    {
        m->free_data();
        for (int i = 0; i < this->size; i++)
        {
            for (int j = 0; j < this->size; j++)
            {
                for (int j = 0; j < this->size; j++)
                {
                   if (fabs(A[i][j]) > 1e-12) m->add(i, j, A[i][j]);
                }
            }
        }
    }

    virtual int get_nnz();

    virtual void print();

    // Return the internal matrix.
    inline scalar **get_A() { return this->A; }

private:
    scalar **A;
};

// **********************************************************************************************************

class CSRMatrix : public Matrix
{
public:
    CSRMatrix(int size);
    CSRMatrix(Matrix *m);
    CSRMatrix(CooMatrix *m);
    CSRMatrix(CSCMatrix *m);
    CSRMatrix(DenseMatrix *m);
    ~CSRMatrix();

    virtual void init();
    virtual void free_data();

    virtual void set_zero()
    {
        _error("CSRMatrix::set_zero() not implemented.");
    }

    void add_from_dense(DenseMatrix *m);
    void add_from_coo(CooMatrix *m);
    void add_from_csc(CSCMatrix *m);

    virtual void add(int m, int n, scalar v)
    {
        _error("CSR matrix add() not implemented.");
    }

    virtual scalar get(int m, int n)
    {
        _error("CSR matrix get() not implemented.");
    }

    virtual int get_size()
    {
        return this->size;
    }
    inline int get_nnz() { return this->nnz; }
    virtual void copy_into(Matrix *m) {
        _error("CSR matrix copy_into() not implemented.");
    }

    virtual void print();

    inline int *get_Ap() { return this->Ap; }
    inline int *get_Ai() { return this->Ai; }
    inline scalar *get_Ax() { return this->Ax; }

private:
    // number of non-zeros
    int nnz;

    int *Ap;
    int *Ai;
    scalar *Ax;
};

// **********************************************************************************************************

class CSCMatrix : public Matrix
{
public:
    CSCMatrix(int size);
    CSCMatrix(Matrix *m);
    CSCMatrix(DenseMatrix *m);
    CSCMatrix(CooMatrix *m);
    CSCMatrix(CSRMatrix *m);
    CSCMatrix(int size, int nnz, int *Ap, int *Ai, scalar *Ax);
    ~CSCMatrix();

    virtual void init();
    virtual void free_data();

    virtual void set_zero()
    {
        _error("CSCMatrix::set_zero() not implemented.");
    }

    void add_from_dense(DenseMatrix *m);
    void add_from_coo(CooMatrix *m);
    void add_from_csr(CSRMatrix *m);

    virtual void add(int m, int n, scalar v)
    {
        _error("CSC matrix add() not implemented.");
    }
    virtual scalar get(int m, int n)
    {
        _error("CSC matrix get() not implemented.");
    }

    virtual int get_size()
    {
        return this->size;
    }
    inline int get_nnz() { return this->nnz;
    }
    virtual void copy_into(Matrix *m)
    {
        _error("CSC matrix copy_into() not implemented.");
    }

    virtual void print();

    inline int *get_Ap() { return this->Ap; }
    inline int *get_Ai() { return this->Ai; }
    inline scalar *get_Ax() { return this->Ax; }

private:
    // number of non-zeros
    int nnz;

    scalar *Ax;

    int *Ap;
    int *Ai;
};

template<typename T>
void dense_to_coo(int size, int nnz, T **Ad, int *row, int *col, T *A);
template<typename T>
void coo_to_csr(int size, int nnz, int *row, int *col, T *A, int *Ap, int *Ai, T *Ax);
template<typename T>
void coo_to_csc(int size, int nnz, int *row, int *col, T *A, int *Ap, int *Ai, T *Ax);
template<typename T>
void csr_to_csc(int size, int nnz, int *Arp, int *Ari, T *Arx, int *Acp, int *Aci, T *Acx);
template<typename T>
void csc_to_csr(int size, int nnz, int *Acp, int *Aci, T *Acx, int *Arp, int *Ari, T *Arx);
template<typename T>
void csc_to_coo(int size, int nnz, int *Ap, int *Ai, T *Ax, int *row, int *col, T *A);
template<typename T>
void csr_to_coo(int size, int nnz, int *Ap, int *Ai, T *Ax, int *row, int *col, T *A);

// matrix vector multiplication
void mat_dot(Matrix *A, scalar *x, scalar *result, int n_dof);
// vector vector multiplication
scalar vec_dot(scalar *r, scalar *s, int n_dof);

void ludcmp(double** a, int n, int* indx, double* d);
void lubksb(double** a, int n, int* indx, double* b);

#endif

