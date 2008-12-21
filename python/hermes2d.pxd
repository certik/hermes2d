cdef extern from "math.h":

    double c_sqrt "sqrt"(double x)

cdef extern from "stdlib.h":
    ctypedef int size_t
    void *malloc (size_t size)

cdef extern from "arrayobject.h":

    ctypedef int intp

    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
        cdef int flags

cdef extern from "Python.h":
    ctypedef void PyObject
    void Py_INCREF(PyObject *x)
    void Py_DECREF(PyObject *x)


cdef extern from "hermes2d.h":

    void hermes2d_initialize(int* argc, char* argv[])
    void hermes2d_finalize()
    int BC_ESSENTIAL "BC_ESSENTIAL"
    int BC_NATURAL "BC_NATURAL"
    int BC_NONE "BC_NONE"
    int c_FN_VAL "FN_VAL"
    int c_FN_DX "FN_DX"
    int c_FN_DY "FN_DY"
    int c_FN_DXX "FN_DXX"
    int c_FN_DYY "FN_DYY"
    int c_FN_DXY "FN_DXY"
    int c_FN_DEFAULT "FN_DEFAULT"
    int c_FN_ALL "FN_ALL"
    int c_EPS_NORMAL "EPS_NORMAL"
    int c_EPS_HIGH "EPS_HIGH"
    int c_verbose_mode "verbose_mode"
    int c_warn_integration "warn_integration"

    cdef struct c_Element "Element":
        int marker
        double get_area()
        double get_diameter()

    cdef struct c_Mesh "Mesh":
        void load(char* filename)
        void save(char* filename)
        void copy(c_Mesh *m)
        void refine_element(int id)
        void refine_all_elements()
        void refine_towards_boundary(int marker, int depth)
        c_Element* get_element(int id)
    c_Mesh *new_Mesh "new Mesh" ()
    void del_Mesh "delete" (c_Mesh *mesh)

    ctypedef struct c_H1Shapeset "H1Shapeset"
    c_H1Shapeset *new_H1Shapeset "new H1Shapeset" ()

    cdef struct c_PrecalcShapeset "PrecalcShapeset"
    c_PrecalcShapeset *new_PrecalcShapeset "new PrecalcShapeset" (c_H1Shapeset *s)

    ctypedef double scalar
    cdef struct EdgePos "EdgePos":
        int marker

    cdef struct c_H1Space "H1Space":
        void set_uniform_order(int tri_order)
        int assign_dofs(int first_dof, int stride)
        void copy_orders(c_H1Space *s, int inc)
        void set_bc_types(int (*bc_type_callback)(int marker))
        void set_bc_values(scalar (*bc_value_callback_by_coord)(int marker,
            double x, double y))
        void set_bc_values_edge "set_bc_values"(scalar (*bc_value_callback_by_edge)(EdgePos *ep))
    c_H1Space *new_H1Space "new H1Space" (c_Mesh *m,
            c_H1Shapeset *h)

    ctypedef struct RealFunction "Function<double>":
        c_Element* get_active_element()
    cdef struct RefMap "RefMap"


    ctypedef scalar (*BiFormFnVol)(RealFunction *fu, RealFunction *fv,
            RefMap *ru, RefMap *rv)

    ctypedef scalar (*LiFormFnVol)(RealFunction *fv, RefMap *rv)

    cdef BiFormFnVol BF_SYM = <BiFormFnVol> 1
    cdef BiFormFnVol BF_ANTISYM = <BiFormFnVol> 2

    cdef struct c_Function "MeshFunction":
        pass

    cdef struct c_MeshFunction "MeshFunction":
        RefMap* get_refmap()

    cdef struct c_Solution "Solution":
        void set_zero(c_Mesh *m)
        void set_fe_solution(c_H1Space *s, c_PrecalcShapeset *pss, scalar *vec)
    c_Solution *new_Solution "new Solution" ()

    cdef struct c_VonMisesFilter "VonMisesFilter"
    c_VonMisesFilter *new_VonMisesFilter "new VonMisesFilter" (c_MeshFunction *sln1, c_MeshFunction *sln2, double l, double m)

    cdef struct c_MagFilter "MagFilter"
    c_MagFilter *new_MagFilter "new MagFilter" (c_MeshFunction *sln1,
            c_MeshFunction *sln2, int item1, int item2)

    cdef struct c_DiffFilter "DiffFilter"
    c_DiffFilter *new_DiffFilter "new DiffFilter" (c_MeshFunction *sln1,
            c_MeshFunction *sln2, int item1, int item2)

    cdef struct c_SumFilter "SumFilter"
    c_SumFilter *new_SumFilter "new SumFilter" (c_MeshFunction *sln1,
            c_MeshFunction *sln2, int item1, int item2)

    cdef struct c_SquareFilter "SquareFilter"
    c_SquareFilter *new_SquareFilter "new SquareFilter" (c_MeshFunction *sln1,
            int item1)


    cdef struct c_DiscreteProblem "DiscreteProblem":
        void set_num_equations(int neq)
        void set_spaces(int n, ...)
        void set_external_fns(int n, ...)
        void set_pss(int n, ...)
        void set_bilinear_form(int i, int j, ...)
        void set_linear_form(int i, ...)
        void create_matrix()
        void set_quiet(int quiet)
        void assemble_matrix_and_rhs()
        void solve_system(int n, ...)
        void save_matrix_matlab(char *filename, char *varname)
        void save_matrix_coo(char *filename)
        void get_matrix(int *Ap, int *Ai, scalar *Ax, int size)
        void copy(c_DiscreteProblem *ep)
        void free_matrix()
    c_DiscreteProblem *new_DiscreteProblem "new DiscreteProblem" ()

    cdef struct c_L2OrthoHP "L2OrthoHP":
        #double calc_error(c_Solution *sln, c_Solution *rsln)
        double calc_error(...)
        void adapt(double thr, int strat)
    c_L2OrthoHP *new_L2OrthoHP "new L2OrthoHP" (int num, ...)

    cdef struct c_H1OrthoHP "H1OrthoHP":
        #double calc_error(c_Solution *sln, c_Solution *rsln)
        double calc_error(...)
        void adapt(double thr, int strat)
    c_H1OrthoHP *new_H1OrthoHP "new H1OrthoHP" (int num, ...)

    double int_u(RealFunction* fu, RefMap* ru)
    double int_l2_norm(RealFunction* fu, RefMap* ru)
    double l2_norm(c_MeshFunction* fu)
    double integrate(c_MeshFunction *sln)
    double int_grad_u_grad_v(RealFunction *fu, RealFunction *fv,
            RefMap *ru, RefMap *rv)
    double int_x_grad_u_grad_v(RealFunction *fu, RealFunction *fv,
            RefMap *ru, RefMap *rv)
    double int_u_dvdx(RealFunction *fu, RealFunction *fv,
            RefMap *ru, RefMap *rv)
    double int_u_dvdy(RealFunction *fu, RealFunction *fv,
            RefMap *ru, RefMap *rv)
    double int_u_v(RealFunction *fu, RealFunction *fv,
            RefMap *ru, RefMap *rv)
    double int_F_u_v(...)
    double int_F_v(...)
    double int_v(RealFunction *fv, RefMap *rv)
    double int_w_nabla_u_v(RealFunction *w1, RealFunction *w2,
            RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv)
    double surf_int_G_v(RealFunction *fv, RefMap *rv, EdgePos *ep)
    double surf_int_v(RealFunction *fv, RefMap *rv, EdgePos *ep)
    double surf_int_u_v(RealFunction *fu, RealFunction *fv,
            RefMap *ru, RefMap *rv, EdgePos *ep)
    double int_a_dudx_dvdx_b_dudy_dvdy(double a, RealFunction *fu, double b,
            RealFunction *fv, RefMap *ru, RefMap *rv)
    double int_a_dudx_dvdy_b_dudy_dvdx(double a, RealFunction *fu, double b,
            RealFunction *fv, RefMap *ru, RefMap *rv)


    cdef struct c_View "View":
        void set_title(char *title)

    cdef struct c_ScalarView "ScalarView":
        void show(c_MeshFunction *s, ...)
        void show_mesh(int show)
        void show_scale(int show)
        void set_min_max_range(double min, double max)
        void set_title(char *title)
    c_ScalarView *new_ScalarView "new ScalarView" (char *title, ...)

    cdef struct c_BaseView "BaseView":
        void show(c_H1Space *s)
    c_BaseView *new_BaseView "new BaseView" (char *title, ...)

    cdef struct c_MeshView "MeshView":
        void show(c_Mesh *s)
    c_MeshView *new_MeshView "new MeshView" (char *title, ...)

    cdef struct c_OrderView "OrderView":
        void show(c_H1Space *s)
    c_OrderView *new_OrderView "new OrderView" (char *title, ...)

    cdef struct c_VectorView "VectorView":
        void show(c_MeshFunction *s1, c_MeshFunction *s2, double eps)
        void show_scale(int show)
        void set_min_max_range(double min, double max)
    c_VectorView *new_VectorView "new VectorView" (char *title, ...)

    cdef struct c_MatrixView "MatrixView":
        void show(c_DiscreteProblem *ep)
    c_MatrixView *new_MatrixView "new MatrixView" (char *title, ...)


cdef class DiscreteProblem:
    cdef c_DiscreteProblem *thisptr

cdef class H1Space:
    cdef c_H1Space *thisptr

cdef class Transformable:
    pass

cdef class Function(Transformable):
    cdef c_Function *thisptr

cdef class ScalarFunction(Function):
    pass

cdef class MeshFunction(ScalarFunction):
    pass

cdef class Solution(MeshFunction):
    pass
