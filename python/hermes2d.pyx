FN_DX = c_FN_DX
FN_DY = c_FN_DY
FN_VAL = c_FN_VAL
FN_DX = c_FN_DX
FN_DY = c_FN_DY
FN_DXX = c_FN_DXX
FN_DYY = c_FN_DYY
FN_DXY = c_FN_DXY
FN_DEFAULT = c_FN_DEFAULT
FN_ALL = c_FN_ALL
EPS_NORMAL = c_EPS_NORMAL
EPS_HIGH = c_EPS_HIGH

cdef class Element:
    cdef c_Element *thisptr

    def get_diameter(self):
        return self.thisptr.get_diameter()

    def get_area(self):
        return self.thisptr.get_area()

cdef class Mesh:
    cdef c_Mesh *thisptr

    def __cinit__(self):
        self.thisptr = new_Mesh()

    def __dealloc__(self):
        del_Mesh(self.thisptr)

    def copy(self, Mesh m):
        self.thisptr.copy(m.thisptr)

    def load(self, char* filename):
        self.thisptr.load(filename)

    def save(self, char* filename):
        self.thisptr.save(filename)

    def refine_element(self, int id):
        self.thisptr.refine_element(id)

    def refine_all_elements(self):
        self.thisptr.refine_all_elements()

    def refine_towards_boundary(self, int marker, int depth):
        self.thisptr.refine_towards_boundary(marker, depth)

    def get_element(self, int id):
        cdef Element e = Element()
        e.thisptr = self.thisptr.get_element(id)
        return e

cdef class H1Shapeset:
    cdef c_H1Shapeset *thisptr

    def __cinit__(self):
        self.thisptr = new_H1Shapeset()

cdef class PrecalcShapeset:
    cdef c_PrecalcShapeset *thisptr

    def __cinit__(self, H1Shapeset h):
        self.thisptr = new_PrecalcShapeset(h.thisptr)

cdef class H1Space:

    def __cinit__(self, Mesh m, H1Shapeset s):
        self.thisptr = new_H1Space(m.thisptr, s.thisptr)

    def set_uniform_order(self, int tri_order):
        self.thisptr.set_uniform_order(tri_order)

    def assign_dofs(self, first_dof=0, stride=1):
        return self.thisptr.assign_dofs(first_dof, stride)

    def copy_orders(self, H1Space s, int inc=0):
        self.thisptr.copy_orders(s.thisptr, inc)

cdef class Transformable:
    pass

cdef class Function(Transformable):
    pass

cdef class ScalarFunction(Function):
    pass

cdef class MeshFunction(ScalarFunction):

    def __add__(x, y):
        return SumFilter(x, y)

    def __sub__(x, y):
        return DiffFilter(x, y)

    def __pow__(x, y, z):
        if y == 2:
            return SquareFilter(x)
        return NotImplemented

    def __neg__(x):
        return x-x-x

    #def integrate(self):
    #    cdef c_MeshFunction *m = <c_MeshFunction *>(self.thisptr)
    #    return integrate(m)

    def l2_norm(self):
        cdef c_MeshFunction *m = <c_MeshFunction *>(self.thisptr)
        return l2_norm(m)

    def get_pt_value(self, x, y):
        cdef c_MeshFunction *m = <c_MeshFunction *>(self.thisptr)
        return m.get_pt_value(x, y)

cdef class Solution(MeshFunction):

    def __cinit__(self):
        self.thisptr = <c_Function *>new_Solution()

    def set_zero(self, Mesh m):
        (<c_Solution *>(self.thisptr)).set_zero(m.thisptr)

    def set_fe_solution(self, H1Space s, PrecalcShapeset pss, ndarray v):
        cdef int n = len(v)
        cdef int i
        from numpy import array
        cdef ndarray vec = array(v, dtype="double")
        cdef scalar *pvec = <scalar *>vec.data
        # XXX: memory leak here?
        cdef scalar *p = <scalar *>malloc(sizeof(scalar)*n)
        for i in range(n):
            p[i] = pvec[i]
        (<c_Solution *>(self.thisptr)).set_fe_solution(s.thisptr, pss.thisptr, p)

cdef class Filter(MeshFunction):
    pass

cdef class SimpleFilter(Filter):
    pass

cdef class VonMisesFilter(Filter):

    def __cinit__(self, MeshFunction sln1, MeshFunction sln2, double l, double m):
        self.thisptr = <c_Function *>new_VonMisesFilter(<c_MeshFunction *>sln1.thisptr, <c_MeshFunction *>sln2.thisptr, l, m)

cdef class MagFilter(SimpleFilter):

    def __cinit__(self, MeshFunction sln1, MeshFunction sln2, int item1, int item2):
        self.thisptr = <c_Function *>new_MagFilter(<c_MeshFunction *>sln1.thisptr, <c_MeshFunction *>sln2.thisptr, item1, item2)

cdef class DiffFilter(SimpleFilter):

    def __cinit__(self, MeshFunction sln1, MeshFunction sln2,
            int item1=FN_VAL, int item2=FN_VAL):
        self.thisptr = <c_Function *>new_DiffFilter(<c_MeshFunction *>sln1.thisptr, <c_MeshFunction *>sln2.thisptr, item1, item2)

cdef class SumFilter(SimpleFilter):

    def __cinit__(self, MeshFunction sln1, MeshFunction sln2,
            int item1=FN_VAL, int item2=FN_VAL):
        self.thisptr = <c_Function *>new_SumFilter(<c_MeshFunction *>sln1.thisptr, <c_MeshFunction *>sln2.thisptr, item1, item2)

cdef class SquareFilter(SimpleFilter):

    def __cinit__(self, MeshFunction sln, int item=FN_VAL):
        self.thisptr = <c_Function *>new_SquareFilter(<c_MeshFunction *>sln.thisptr, item)

cdef class DiscreteProblem:

    def __cinit__(self):
        self.thisptr = new_DiscreteProblem()

    def copy(self, DiscreteProblem ep):
        self.thisptr.copy(ep.thisptr)

    def set_num_equations(self, int neq):
        self.thisptr.set_num_equations(neq)

    def set_spaces(self, *args):
        cdef int n = len(args)
        cdef H1Space a, b, c
        if n == 1:
            a = args[0]
            self.thisptr.set_spaces(n, a.thisptr)
        elif n == 2:
            a, b = args
            self.thisptr.set_spaces(n, a.thisptr, b.thisptr)
        elif n == 3:
            a, b, c = args
            self.thisptr.set_spaces(n, a.thisptr, b.thisptr, c.thisptr)
        else:
            raise NotImplementedError()

    def set_external_fns(self, *args):
        cdef int n = len(args)
        cdef Solution a, b, c
        if n == 1:
            a = args[0]
            self.thisptr.set_external_fns(n, a.thisptr)
        elif n == 2:
            a, b = args
            self.thisptr.set_external_fns(n, a.thisptr, b.thisptr)
        elif n == 3:
            a, b, c = args
            self.thisptr.set_external_fns(n, a.thisptr, b.thisptr, c.thisptr)
        else:
            raise NotImplementedError()

    def set_pss(self, *args):
        cdef int n = len(args)
        cdef PrecalcShapeset s
        if n == 1:
            s = args[0]
            self.thisptr.set_pss(n, s.thisptr)
        else:
            raise NotImplementedError()

    #def set_bilinear_form(self, int i, int j, BiFormFnVol unsym):
    #    self.thisptr.set_bilinear_form(i, j, unsym)

    def create_matrix(self):
        self.thisptr.create_matrix()

    def set_quiet(self, quiet=True):
        self.thisptr.set_quiet(quiet)

    def assemble_matrix_and_rhs(self):
        self.thisptr.assemble_matrix_and_rhs()

    def solve_system(self, *args):
        cdef int n = len(args)
        cdef Solution a, b, c
        if n == 1:
            a = args[0]
            self.thisptr.solve_system(n, a.thisptr)
        elif n == 2:
            a, b = args
            self.thisptr.solve_system(n, a.thisptr, b.thisptr)
        elif n == 3:
            a, b, c = args
            self.thisptr.solve_system(n, a.thisptr, b.thisptr, c.thisptr)
        else:
            raise NotImplementedError()

    def save_matrix_matlab(self, char *filename, char *varname):
        self.thisptr.save_matrix_matlab(filename, varname)

    def save_matrix_coo(self, char *filename):
        self.thisptr.save_matrix_coo(filename)

    def get_matrix_csc(self):
        """
        Returns the matrix A as a (Ap, Ai, Ax) tuple.

        See also DiscreteProblem.get_matrix() to get a SciPy matrix.
        """
        cdef int *Ap, *Ai, n, nnz
        cdef scalar *Ax
        self.thisptr.get_matrix(Ap, Ai, Ax, n)
        nnz = Ap[n]
        from numpy import zeros
        cdef ndarray aAp = zeros(n+1, dtype="int32")
        cdef ndarray aAi = zeros(nnz, dtype="int32")
        cdef ndarray aAx = zeros(nnz, dtype="double")
        cdef int *pAp = <int *>aAp.data
        cdef int *pAi = <int *>aAi.data
        cdef double *pAx = <double *>aAx.data
        cdef int i
        for i in range(n+1):
            pAp[i] = Ap[i]
        for i in range(nnz):
            pAi[i] = Ai[i]
            pAx[i] = Ax[i]
        return aAp, aAi, aAx

    def get_matrix(self):
        """
        Returns the global matrix A as a SciPy matrix.
        """
        from scipy.sparse import csc_matrix
        Ap, Ai, Ax = self.get_matrix_csc()
        return csc_matrix((Ax, Ai, Ap))

    def free_matrix(self):
        self.thisptr.free_matrix()

cdef class L2OrthoHP:
    cdef c_L2OrthoHP *thisptr

    def __cinit__(self, *args):
        cdef int n = len(args)
        cdef H1Space a, b, c
        if n == 1:
            a = args[0]
            self.thisptr = new_L2OrthoHP(n, a.thisptr)
        elif n == 2:
            a, b = args
            self.thisptr = new_L2OrthoHP(n, a.thisptr, b.thisptr)
        elif n == 3:
            a, b, c = args
            self.thisptr = new_L2OrthoHP(n, a.thisptr, b.thisptr, c.thisptr)
        else:
            raise NotImplementedError()

    def calc_error(self, MeshFunction sln, MeshFunction rsln):
        return self.thisptr.calc_error(sln.thisptr, rsln.thisptr)

    def adapt(self, double thr, int strat=0):
        self.thisptr.adapt(thr, strat)

cdef class H1OrthoHP:
    cdef c_H1OrthoHP *thisptr

    def __cinit__(self, *args):
        cdef int n = len(args)
        cdef H1Space a, b, c
        if n == 1:
            a = args[0]
            self.thisptr = new_H1OrthoHP(n, a.thisptr)
        elif n == 2:
            a, b = args
            self.thisptr = new_H1OrthoHP(n, a.thisptr, b.thisptr)
        elif n == 3:
            a, b, c = args
            self.thisptr = new_H1OrthoHP(n, a.thisptr, b.thisptr, c.thisptr)
        else:
            raise NotImplementedError()

    def calc_error(self, MeshFunction sln, MeshFunction rsln):
        return self.thisptr.calc_error(sln.thisptr, rsln.thisptr)

    def adapt(self, double thr, int strat=0):
        self.thisptr.adapt(thr, strat)

cdef class View:
    pass

cdef class ScalarView(View):
    cdef c_ScalarView *thisptr

    def __cinit__(self, char *title, int x=-1, int y=-1, int width=1000,
                int height = 800):
        self.thisptr = new_ScalarView(title, x, y, width, height)

    def show(self, MeshFunction sln, eps=None):
        cdef Function s = <Function>sln
        if eps is None:
            self.thisptr.show(<c_MeshFunction *>s.thisptr)
        else:
            self.thisptr.show(<c_MeshFunction *>s.thisptr, <double>eps)

    def show_mesh(self, show=True):
        self.thisptr.show_mesh(show)

    def show_scale(self, show=True):
        self.thisptr.show_scale(show)

    def set_min_max_range(self, double min, double max):
        self.thisptr.set_min_max_range(min, max)

    def set_title(self, char *title):
        ((self.thisptr)).set_title(title)
        #(<c_View *>(self.thisptr)).set_title(title)

cdef class BaseView(View):
    cdef c_BaseView *thisptr

    def __cinit__(self, char *title="BaseView", int x=-1, int y=-1,
            int width=1000, int height = 800):
        self.thisptr = new_BaseView(title, x, y, width, height)

    def show(self, H1Space s):
        self.thisptr.show(s.thisptr)

cdef class MeshView(View):
    cdef c_MeshView *thisptr

    def __cinit__(self, char *title="MeshView", int x=-1, int y=-1,
            int width=1000, int height = 800):
        self.thisptr = new_MeshView(title, x, y, width, height)

    def show(self, Mesh s):
        self.thisptr.show(s.thisptr)

cdef class OrderView(View):
    cdef c_OrderView *thisptr

    def __cinit__(self, char *title="OrderView", int x=-1, int y=-1,
            int width=1000, int height = 800):
        self.thisptr = new_OrderView(title, x, y, width, height)

    def show(self, H1Space s):
        self.thisptr.show(s.thisptr)

cdef class VectorView(View):
    cdef c_VectorView *thisptr

    def __cinit__(self, char *title, int x=-1, int y=-1, int width=1000,
            int height = 800):
        self.thisptr = new_VectorView(title, x, y, width, height)

    def show(self, MeshFunction s1, MeshFunction s2, eps = 0.008):
        (<c_VectorView *>(self.thisptr)).show(<c_MeshFunction *>s1.thisptr, <c_MeshFunction *>s2.thisptr, eps)

    def show_scale(self, show=True):
        self.thisptr.show_scale(show)

    def set_min_max_range(self, double min, double max):
        self.thisptr.set_min_max_range(min, max)

cdef class MatrixView(View):
    cdef c_MatrixView *thisptr

    def __cinit__(self, char *title, int x=-1, int y=-1, int width=1000,
            int height = 800):
        self.thisptr = new_MatrixView(title, x, y, width, height)

    def show(self, DiscreteProblem ep):
        self.thisptr.show(ep.thisptr)

def initialize():
    cdef int argc = 1
    cdef char* argv = "poisson"
    hermes2d_initialize(&argc, &argv)

def set_verbose(verbose):
    """
    Sets the global verbose_mode variable.

    That variable controls how verbose hermes2d is.

    Returns the old status.
    """
    global c_verbose_mode
    old_flag = c_verbose_mode
    c_verbose_mode = verbose
    return old_flag

def set_warn_integration(warn_integration):
    """
    Sets the global warn_integration variable.

    That variable controls if warn_order() in the C++ hermes2d should emit the
    warning about integration rules.

    Returns the old status.
    """
    global c_warn_integration
    old_flag = c_warn_integration
    c_warn_integration = warn_integration
    return old_flag

def finalize():
    hermes2d_finalize()
