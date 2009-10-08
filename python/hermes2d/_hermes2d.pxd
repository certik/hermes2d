cdef extern from "math.h":

    double c_sqr "sqr"(double x)
    double c_sqrt "sqrt"(double x)
    double c_atan "atan"(double x)
    double c_pi "M_PI"

cdef extern from "stdlib.h":

    ctypedef unsigned long size_t
    void *malloc (size_t size)
    void free(void *mem)
    void *memcpy(void *dst, void *src, long n)

    void exit(int exit_code)

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

cdef extern from "stdcython.h":
    void init_global_empty_tuple()
    object PY_NEW(object t)

cdef extern from "hermes2d.h":

    # This is just the C++ "delete" statement
    void delete(...)

    void hermes2d_initialize(int* argc, char* argv[])
    void hermes2d_finalize()
    void throw_exception(char *msg)
    #void finish_glut_main_loop(int force_quit)
    int BC_ESSENTIAL "BC_ESSENTIAL"
    int BC_NATURAL "BC_NATURAL"
    int BC_NONE "BC_NONE"
    int ANTISYM "ANTISYM"
    int UNSYM "UNSYM"
    int SYM "SYM"
    int ANY "ANY"
    int c_FN_VAL "FN_VAL"
    int c_FN_VAL_0 "FN_VAL_0"
    int c_FN_DX "FN_DX"
    int c_FN_DY "FN_DY"
    int c_FN_DXX "FN_DXX"
    int c_FN_DYY "FN_DYY"
    int c_FN_DXY "FN_DXY"
    int c_FN_DEFAULT "FN_DEFAULT"
    int c_FN_ALL "FN_ALL"
    int c_EPS_LOW "EPS_LOW"
    int c_EPS_NORMAL "EPS_NORMAL"
    int c_EPS_HIGH "EPS_HIGH"
    int c_verbose_mode "verbose_mode"
    int c_info_mode "info_mode"
    int c_warn_integration "warn_integration"

    ctypedef double double4[4]
    ctypedef double double3[3]
    ctypedef int int3[3]
    ctypedef int int2[2]

    cdef struct c_CurvMap "CurvMap":
        int toplevel
        int order

    cdef struct c_Node "Node":
        int id
        unsigned ref
        unsigned type
        unsigned bnd
        unsigned used
        double x, y

    cdef struct c_Element "Element":
        int id
        unsigned nvert
        unsigned active
        unsigned used
        int marker
        int userdata
        int iro_cache
        c_Node* vn[4]
        c_Node* en[4]
        c_Element* sons[4]
        c_CurvMap* cm
        double get_area()
        double get_diameter()

    cdef struct c_Mesh "Mesh":
        void load(char* filename)
        void load_str(char* mesh)
        void save(char* filename)
        void copy(c_Mesh *m)
        void refine_element(int id)
        void refine_all_elements()
        void refine_towards_boundary(int marker, int depth)
        c_Element* get_element(int id)
        int get_num_elements()
        int get_num_base_elements()
        int get_num_active_elements()
        int get_max_element_id()

    c_Mesh *new_Mesh "new Mesh" ()

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
        int get_element_order(int id)
        void set_bc_types(int (*bc_type_callback)(int marker))
        void set_bc_values(scalar (*bc_value_callback_by_coord)(int marker,
            double x, double y))
        void set_bc_values_edge "set_bc_values"(scalar (*bc_value_callback_by_edge)(EdgePos *ep))
    c_H1Space *new_H1Space "new H1Space" (c_Mesh *m,
            c_H1Shapeset *h)

    ctypedef struct RealFunction "Function<double>":
        c_Element* get_active_element()
    cdef struct RefMap "RefMap"
    ctypedef struct c_ScalarFunction "Function<scalar>"


    ctypedef scalar (*BiFormFnVol)(RealFunction *fu, RealFunction *fv,
            RefMap *ru, RefMap *rv)

    ctypedef scalar (*LiFormFnVol)(RealFunction *fv, RefMap *rv)

    cdef BiFormFnVol BF_SYM = <BiFormFnVol> 1
    cdef BiFormFnVol BF_ANTISYM = <BiFormFnVol> 2

    cdef struct c_Function "MeshFunction":
        pass

    cdef struct c_MeshFunction "MeshFunction":
        RefMap* get_refmap()
        c_Mesh* get_mesh()
        #scalar get_pt_value(double x, double y, int item)
        scalar get_pt_value(double x, double y)

    cdef struct c_Solution "Solution":
        void set_zero(c_Mesh *m)
        void set_fe_solution(c_H1Space *s, c_PrecalcShapeset *pss, scalar *vec)
        void get_fe_solution(int *Ylen, scalar **Y)
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

    cdef struct c_WeakForm "WeakForm":
        void add_biform(int i, int j, ...)
        void add_biform_surf(int i, int j, ...)
        void add_liform(int i, ...)
        void add_liform_data(int i, void *data)
        void add_liform_surf(int i, ...)
    c_WeakForm *new_WeakForm "new WeakForm" (int neq)

    cdef struct c_Solver "Solver":
        pass

    cdef struct c_LinSystem "LinSystem":
        void set_spaces(int n, ...)
        void set_pss(int n, ...)
        c_H1Space *get_space(int n)
        c_PrecalcShapeset *get_pss(int n)
        void copy(c_LinSystem *sys)
        void assemble()
        int solve(int n, ...)
        void save_matrix_matlab(char *filename, char *varname)
        void get_matrix(int *Ap, int *Ai, scalar *Ax, int size)
        void get_rhs(scalar *RHS, int size)
    c_LinSystem *new_LinSystem "new LinSystem" (c_WeakForm *wf,
            c_Solver *solver)

    cdef struct c_RefSystem "RefSystem":
        void assemble()
        c_H1Space *get_ref_space(int eq)
    c_RefSystem *new_RefSystem "new RefSystem" (c_LinSystem *ls)

    #cdef struct c_DiscreteProblem "DiscreteProblem":
    #    void set_num_equations(int neq)
    #    void set_external_fns(int n, ...)
    #    void set_bilinear_form(int i, int j, ...)
    #    void set_linear_form(int i, ...)
    #    void create_matrix()
    #    void set_quiet(int quiet)
    #    void assemble_matrix_and_rhs()
    #    void solve_system(int n, ...)
    #    void save_matrix_coo(char *filename)
    #    void free_matrix()
    #c_DiscreteProblem *new_DiscreteProblem "new DiscreteProblem" ()

    cdef struct c_L2OrthoHP "L2OrthoHP":
        #double calc_error(c_Solution *sln, c_Solution *rsln)
        double calc_error(...)
        double calc_error_n(int n, ...)
        void adapt(double thr, int strat, int h_only)
    c_L2OrthoHP *new_L2OrthoHP "new L2OrthoHP" (int num, ...)

    cdef struct c_H1OrthoHP "H1OrthoHP":
        #double calc_error(c_Solution *sln, c_Solution *rsln)
        int num
        double calc_error(...)
        double calc_error_n(int n, ...)
        void adapt(double thr, int strat, int h_only)
    c_H1OrthoHP *new_H1OrthoHP "new H1OrthoHP" (int num, ...)

    cdef struct c_Linearizer "Linearizer":
        void process_solution(c_MeshFunction* sln, ...)
        void lock_data()
        void unlock_data()
        double3* get_vertices()
        int get_num_vertices()
        int3* get_triangles()
        int get_num_triangles()
        int3* get_edges()
        int get_num_edges()
        double get_min_value()
        double get_max_value()
        void save_data(char* filename)
        void load_data(char* filename)
    c_Linearizer *new_Linearizer "new Linearizer" ()

    cdef struct c_Vectorizer "Vectorizer":
        void process_solution(c_MeshFunction* xsln, ...)
        double4* get_vertices()
        int2* get_dashes()
        int get_num_dashes()
    c_Vectorizer *new_Vectorizer "new Vectorizer" ()


    double int_u(RealFunction* fu, RefMap* ru)
    double int_l2_norm(RealFunction* fu, RefMap* ru)
    double l2_norm(c_MeshFunction* fu)
    double h1_norm(c_MeshFunction* fu)
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
    void View_wait "View::wait"()

    cdef struct c_ScalarView "ScalarView":
        void show(c_MeshFunction *s, ...)
        void show_mesh(int show)
        void show_scale(int show)
        void set_min_max_range(double min, double max)
        void set_title(char *title)
        void wait()
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

    #cdef struct c_MatrixView "MatrixView":
    #    pass
        #void show(c_DiscreteProblem *ep)
    #c_MatrixView *new_MatrixView "new MatrixView" (char *title, ...)

cdef extern from "dummy_solver.h":

    cdef struct c_DummySolver "DummySolver":
        pass
    c_DummySolver *new_DummySolver "new DummySolver" ()


cdef class LinSystem:
    cdef c_LinSystem *thisptr
    cdef object _spaces
    cdef object _pss

cdef class RefSystem(LinSystem):
    pass

cdef class Solver:
    cdef c_Solver *thisptr

cdef class DummySolver(Solver):
    pass

cdef class WeakForm:
    cdef c_WeakForm *thisptr

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

cdef class Linearizer:
    cdef c_Linearizer *thisptr

cdef class Vectorizer(Linearizer):
    pass
