# Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
# Distributed under the terms of the BSD license (see the LICENSE
# file for the exact terms).
# Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

cdef inline object cp2str(const_char_p p):
    if p == NULL: return None
    else:         return p

cdef inline char_p str2cp(object s) except ? NULL:
    if s is None: return NULL
    else:         return s


#-----------------------------------------------------------------------
# Matrix classes:

cdef class Matrix:
    cdef c_Matrix *thisptr

    def get_size(self):
        return self.thisptr.get_size()

    def add(self, int m, int n, double v):
        self.thisptr.add(m, n, v)

cdef class SparseMatrix(Matrix):
    pass

cdef class CooMatrix(SparseMatrix):

    def __init__(self, size=0, is_complex=False):
        self.thisptr = <c_Matrix *>new_CooMatrix(size, is_complex)

    def add(self, int m, int n, v):
        if self.thisptr.is_complex():
            self.thisptr.add_cplx(m, n, v)
        else:
            self.thisptr.add(m, n, v)

    @property
    def row_col_data(self):
        """
        Returns (row, col, data) arrays.
        """
        from numpy import empty
        cdef c_CooMatrix *_thisptr = <c_CooMatrix*>(self.thisptr)
        cdef int n, len
        cdef int *crow, *ccol
        cdef double *cdata
        cdef cplx *ccdata
        if self.thisptr.is_complex():
            len = _thisptr.triplets_len_cplx()
        else:
            len = _thisptr.triplets_len()
        row = empty([len], dtype="int32")
        numpy2c_int_inplace(row, &crow, &n)
        col = empty([len], dtype="int32")
        numpy2c_int_inplace(col, &ccol, &n)
        if self.thisptr.is_complex():
            data = empty([len], dtype="complex128")
            numpy2c_double_complex_inplace(data, &ccdata, &n)
            _thisptr.get_row_col_data_cplx(crow, ccol, ccdata)
            return row, col, data
        else:
            data = empty([len], dtype="double")
            numpy2c_double_inplace(data, &cdata, &n)
            _thisptr.get_row_col_data(crow, ccol, cdata)
        return row, col, data

    def to_scipy_coo(self):
        """
        Converts itself to the scipy sparse COO format.
        """
        from scipy.sparse import coo_matrix
        row, col, data = self.row_col_data
        n = self.get_size()
        return coo_matrix((data, (row, col)), shape=(n, n))

    def __str__(self):
        return str(self.to_scipy_coo())

cdef class CSRMatrix(SparseMatrix):

    def __init__(self, M):
        if isinstance(M, (int, long)):
            size = M
            self.thisptr = <c_Matrix *>new_CSRMatrix_size(size)
        elif isinstance(M, CooMatrix):
            self.thisptr = <c_Matrix *>new_CSRMatrix_coo_matrix(
                    <c_CooMatrix*>(py2c_Matrix(M).thisptr))
        elif isinstance(M, CSCMatrix):
            self.thisptr = <c_Matrix *>new_CSRMatrix_csc_matrix(
                    <c_CSCMatrix*>(py2c_Matrix(M).thisptr))
        else:
            raise Exception("Not implemented.")

    @property
    def IA(self):
        """
        Returns (row, col, data) arrays.
        """
        cdef c_CSRMatrix *_thisptr = <c_CSRMatrix*>(self.thisptr)
        return c2numpy_int_inplace(_thisptr.get_IA(), self.get_size()+1)

    @property
    def JA(self):
        """
        Returns (row, col, data) arrays.
        """
        cdef c_CSRMatrix *_thisptr = <c_CSRMatrix*>(self.thisptr)
        return c2numpy_int_inplace(_thisptr.get_JA(), _thisptr.get_nnz())

    @property
    def A(self):
        """
        Returns (row, col, data) arrays.
        """
        cdef c_CSRMatrix *_thisptr = <c_CSRMatrix*>(self.thisptr)
        if self.thisptr.is_complex():
            return c2numpy_double_complex_inplace(_thisptr.get_A_cplx(),
                    _thisptr.get_nnz())
        else:
            return c2numpy_double_inplace(_thisptr.get_A(), _thisptr.get_nnz())

    def to_scipy_csr(self):
        """
        Converts itself to the scipy sparse CSR format.
        """
        from scipy.sparse import csr_matrix
        n = self.get_size()
        return csr_matrix((self.A, self.JA, self.IA), shape=(n, n))

    def __str__(self):
        return str(self.to_scipy_csr())

cdef class CSCMatrix(SparseMatrix):

    def __init__(self, M):
        if isinstance(M, (int, long)):
            size = M
            self.thisptr = <c_Matrix *>new_CSCMatrix_size(size)
        elif isinstance(M, CooMatrix):
            self.thisptr = <c_Matrix *>new_CSCMatrix_coo_matrix(
                    <c_CooMatrix*>(py2c_Matrix(M).thisptr))
        elif isinstance(M, CSRMatrix):
            self.thisptr = <c_Matrix *>new_CSCMatrix_csr_matrix(
                    <c_CSRMatrix*>(py2c_Matrix(M).thisptr))
        else:
            raise Exception("Not implemented.")

    @property
    def IA(self):
        """
        Returns (row, col, data) arrays.
        """
        cdef c_CSCMatrix *_thisptr = <c_CSCMatrix*>(self.thisptr)
        return c2numpy_int_inplace(_thisptr.get_IA(), _thisptr.get_nnz())

    @property
    def JA(self):
        """
        Returns (row, col, data) arrays.
        """
        cdef c_CSCMatrix *_thisptr = <c_CSCMatrix*>(self.thisptr)
        return c2numpy_int_inplace(_thisptr.get_JA(), self.get_size()+1)

    @property
    def A(self):
        """
        Returns (row, col, data) arrays.
        """
        cdef c_CSCMatrix *_thisptr = <c_CSCMatrix*>(self.thisptr)
        if self.thisptr.is_complex():
            return c2numpy_double_complex_inplace(_thisptr.get_A_cplx(),
                    _thisptr.get_nnz())
        else:
            return c2numpy_double_inplace(_thisptr.get_A(), _thisptr.get_nnz())

    def to_scipy_csc(self):
        """
        Converts itself to the scipy sparse CSC format.
        """
        from scipy.sparse import csc_matrix
        n = self.get_size()
        return csc_matrix((self.A, self.IA, self.JA), shape=(n, n))

    def __str__(self):
        return str(self.to_scipy_csc())

cdef Matrix py2c_Matrix(object M):
    return M

cdef api object c2py_CooMatrix(c_CooMatrix *m):
    cdef CooMatrix c
    c = <CooMatrix>PY_NEW(CooMatrix)
    c.thisptr = <c_Matrix *>m
    return c

cdef api object c2py_CSRMatrix(c_CSRMatrix *m):
    cdef CSRMatrix c
    c = <CSRMatrix>PY_NEW(CSRMatrix)
    c.thisptr = <c_Matrix *>m
    return c

cdef api object c2py_CSCMatrix(c_CSCMatrix *m):
    cdef CSCMatrix c
    c = <CSCMatrix>PY_NEW(CSCMatrix)
    c.thisptr = <c_Matrix *>m
    return c

#-----------------------------------------------------------------------
# Common C++ <-> Python+NumPy conversion tools:

import sys
import traceback
# this is important to be called here, otherwise we can't use the NumPy C/API:
import_array()

cdef api object namespace_create():
    return {"verbose": False}

cdef api void namespace_push(object namespace, const_char_p name, object o):
    namespace.update({name: o})

cdef api void namespace_print(object namespace):
    print "-"*80
    print "namespace:"
    print namespace

cdef api object namespace_pull(object namespace, const_char_p name):
    return namespace.get(name)

global_namespace = namespace_create()

cdef api void cmd(const_char_p text):
    """
    Runs the command "text" in the Python namespace.
    """
    run_cmd(text, global_namespace)

cdef api void set_verbose_cmd(int verbose):
    global_namespace["verbose"] = verbose

cdef api void insert_object(const_char_p name, object o):
    """
    Inserts an object into the global namespace.

    Example 1:

    insert_object("a", c2py_int(3));
    cmd("print a");

    This prints "3".

    Example 2:

    int a[3] = {1, 5, 3};
    insert_object("A", c2numpy_int(a, 3));
    cmd("print A");

    This prints "[1  5  3]" (this is how the NumPy array is printed).

    Example 3:

    double a[3] = {1, 5, 3};
    insert_object("A", c2numpy_double(a, 3));
    cmd("print A");

    This prints "[ 1.  5.  3.]" (this is how the NumPy array is printed).
    """
    namespace_push(global_namespace, name, o)

cdef api object get_object(const_char_p name):
    """
    Retrieves an object from the Python namespace.

    Example:

    // put into python:
    int a[3] = {1, 5, 3};
    insert_object("A", c2numpy_int(a, 3));

    // retrieve from python:
    double *A;
    int n;
    numpy2c_double_inplace(get_object("A"), &A, &n);
    """
    return namespace_pull(global_namespace, name)

cdef api object c2py_int(int i):
    return i

cdef api int py2c_int(object i):
    return i

cdef api char* py2c_str(object s):
    return s

cdef api double py2c_double(object i):
    return i

cdef api object c2numpy_int(int *A, int len):
    """
    Construct the integer NumPy array by copying the data.
    """
    from numpy import empty
    cdef ndarray vec = empty([len], dtype="int32")
    cdef int *pvec = <int *>vec.data
    memcpy(pvec, A, len*sizeof(int))
    return vec

cdef api object c2numpy_int_inplace(int *A, int len):
    """
    Construct the integer NumPy array inplace (don't copy any data).
    """
    cdef npy_intp dim = len
    return PyArray_SimpleNewFromData(1, &dim, NPY_INT, A)

cdef api object c2numpy_double(double *A, int len):
    """
    Construct the double NumPy array by copying the data.
    """
    from numpy import empty
    cdef ndarray vec = empty([len], dtype="double")
    cdef double *pvec = <double *>vec.data
    memcpy(pvec, A, len*sizeof(double))
    return vec

cdef api object c2numpy_double_inplace(double *A, int len):
    """
    Construct the double NumPy array inplace (don't copy any data).
    """
    cdef npy_intp dim = len
    return PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, A)

cdef api object c2numpy_double_complex_inplace(cplx *A, int len):
    """
    Construct the double NumPy array inplace (don't copy any data).
    """
    cdef npy_intp dim = len
    return PyArray_SimpleNewFromData(1, &dim, NPY_COMPLEX128, A)

_AA = None

cdef api void numpy2c_int_inplace(object A_n, int **A_c, int *n):
    """
    Returns the C array, that points to the numpy array (inplace).

    Only if strides != sizeof(int), the data get copied first.

    Important note: you need to use the A_c array immediately after calling
    this function in C, otherwise numpy could deallocate the array, especially
    if the _AA global variable was deallocated.
    """
    cdef ndarray A = A_n
    if not (A.nd == 1 and A.strides[0] == sizeof(int)):
        from numpy import array
        A = array(A.flat, dtype="int32")
        # this is needed so that numpy doesn't dealocate the arrays
        global _AA
        _AA = A
    n[0] = len(A)
    A_c[0] = <int *>(A.data)

cdef api void numpy2c_double_inplace(object A_n, double **A_c, int *n):
    """
    Returns the C array, that points to the numpy array (inplace).

    Only if strides != sizeof(double), the data get copied first.

    Important note: you need to use the A_c array immediately after calling
    this function in C, otherwise numpy could deallocate the array, especially
    if the _AA global variable was deallocated.
    """
    cdef ndarray A = A_n
    if not (A.nd == 1 and A.strides[0] == sizeof(double)):
        from numpy import array
        A = array(A.flat, dtype="double")
        # this is needed so that numpy doesn't dealocate the arrays
        global _AA
        _AA = A
    n[0] = len(A)
    A_c[0] = <double *>(A.data)

cdef api void numpy2c_double_complex_inplace(object A_n, cplx **A_c, int *n):
    """
    Returns the C array, that points to the numpy array (inplace).

    Only if strides != sizeof(double), the data get copied first.

    Important note: you need to use the A_c array immediately after calling
    this function in C, otherwise numpy could deallocate the array, especially
    if the _AA global variable was deallocated.
    """
    cdef ndarray A = A_n
    if not (A.nd == 1 and A.strides[0] == sizeof(cplx)):
        from numpy import array
        A = array(A.flat, dtype="complex128")
        # this is needed so that numpy doesn't dealocate the arrays
        global _AA
        _AA = A
    n[0] = len(A)
    A_c[0] = <cplx *>(A.data)

cdef api void run_cmd(const_char_p text, object namespace):
    try:
        verbose = namespace.get("verbose")
        if verbose:
            print "got a text:", text
        if verbose:
            print "evaluting in the namespace:"
            print namespace
        code = compile(text, "", "exec")
        eval(code, {}, namespace)
        if verbose:
            print "new namespace:"
            print namespace
    except SystemExit, e:
        try:
            exit_code = int(e)
        except:
            exit_code = -1
        exit(exit_code)
    except:
        etype, value, tb = sys.exc_info()
        s = "".join(traceback.format_exception(etype, value, tb))
        s = "Exception raised in the Python code:\n" + s
        throw_exception(s)

def init_hermes2d_wrappers():
    init_global_empty_tuple()

init_global_empty_tuple()
