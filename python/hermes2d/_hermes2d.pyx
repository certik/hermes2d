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

cdef class Node:
    cdef c_Node *thisptr

    @property
    def id(self):
        return self.thisptr.id

    @property
    def ref(self):
        return self.thisptr.ref

    @property
    def type(self):
        return self.thisptr.type

    @property
    def bnd(self):
        return self.thisptr.bnd

    @property
    def used(self):
        return self.thisptr.used == 1

    def __str__(self):
        return "Node %d: coord=%r, used=%r" % (self.id, self.coord, self.used)

    @property
    def coord(self):
        return self.thisptr.x, self.thisptr.y

cdef class Element:
    cdef c_Element *thisptr

    @property
    def id(self):
        return self.thisptr.id

    @property
    def nvert(self):
        return self.thisptr.nvert

    @property
    def active(self):
        return self.thisptr.active == 1

    @property
    def used(self):
        return self.thisptr.used == 1

    @property
    def marker(self):
        return self.thisptr.marker

    @property
    def nodes_vertex(self):
        return [self.get_vertex_node(i) for i in range(self.nvert)]

    @property
    def nodes_edge(self):
        return [self.get_edge_node(i) for i in range(self.nvert)]

    def get_vertex_node(self, int id):
        cdef Node n = Node()
        n.thisptr = self.thisptr.vn[id]
        return n

    def get_edge_node(self, int id):
        cdef Node n = Node()
        n.thisptr = self.thisptr.en[id]
        return n

    def get_son_element(self, int id):
        cdef Element e = Element()
        e.thisptr = self.thisptr.sons[id]
        return e

    def __str__(self):
        nodes_id = [node.id for node in self.nodes_vertex]
        return "Element %d: nodes_id=%r, active=%r, marker=%d, used=%r" % \
                (self.id, nodes_id, self.active, self.marker, self.used)

    def get_diameter(self):
        return self.thisptr.get_diameter()

    def get_area(self):
        return self.thisptr.get_area()

cdef class Mesh:
    cdef c_Mesh *thisptr

    def __cinit__(self):
        self.thisptr = new_Mesh()

    def __dealloc__(self):
        delete(self.thisptr)

    def copy(self, Mesh m):
        self.thisptr.copy(m.thisptr)

    def load(self, char* filename):
        self.thisptr.load(filename)

    def load_str(self, char* mesh):
        self.thisptr.load_str(mesh)

    @property
    def num_elements(self):
        return self.thisptr.get_num_elements()

    @property
    def num_base_elements(self):
        return self.thisptr.get_num_base_elements()

    @property
    def num_active_elements(self):
        return self.thisptr.get_num_active_elements()

    @property
    def max_element_id(self):
        return self.thisptr.get_max_element_id()

    def __str__(self):
        return "Mesh with %d elements: active=%d, base=%d, max_id=%d" % \
                (self.num_elements, self.num_active_elements,
                        self.num_base_elements, self.max_element_id)

    def create(self, nodes, elements, boundary, nurbs):
        """
        Creates a mesh from a list of nodes, elements, boundary and nurbs.

        Example:

        >>> m = hermes2d.Mesh()
        >>> m.create([
                    [0, -1],
                    [1, -1],
                    [-1, 0],
                    [0, 0],
                    [1, 0],
                    [-1, 1],
                    [0, 1],
                    [0.707106781, 0.707106781],
                ], [
                    [0, 1, 4, 3, 0],
                    [3, 4, 7, 0],
                    [3, 7, 6, 0],
                    [2, 3, 6, 5, 0],
                ], [
                    [0, 1, 1],
                    [1, 4, 2],
                    [3, 0, 4],
                    [4, 7, 2],
                    [7, 6, 2],
                    [2, 3, 4],
                    [6, 5, 2],
                    [5, 2, 3],
                ], [
                    [4, 7, 45],
                    [7, 6, 45],
                ])

        """
        m = "1 0\n"
        m += "%d\n" % len(nodes)
        for node in nodes:
            m += "%f %f\n" % tuple(node)
        m += "%d\n" % len(elements)
        for el in elements:
            m += ("%d "*len(el)) % tuple(el)
            m += "\n"
        m += "%d\n" % len(boundary)
        for b in boundary:
            m += "%d %d %d\n" % tuple(b)
        m += "%d\n" % len(nurbs)
        for n in nurbs:
            m += "%d %d\n0\n%f\n" % tuple(n)
        m += "0\n"
        self.load_str(m)

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

    def __dealloc__(self):
        delete(self.thisptr)

cdef class PrecalcShapeset:
    cdef c_PrecalcShapeset *thisptr

    def __cinit__(self, H1Shapeset h):
        self.thisptr = new_PrecalcShapeset(h.thisptr)

    def __dealloc__(self):
        delete(self.thisptr)

cdef class H1Space:

    def __init__(self, Mesh m, H1Shapeset s):
        self.thisptr = new_H1Space(m.thisptr, s.thisptr)

    def __dealloc__(self):
        delete(self.thisptr)

    def set_uniform_order(self, int tri_order):
        self.thisptr.set_uniform_order(tri_order)

    def assign_dofs(self, first_dof=0, stride=1):
        return self.thisptr.assign_dofs(first_dof, stride)

    def copy_orders(self, H1Space s, int inc=0):
        self.thisptr.copy_orders(s.thisptr, inc)

cdef H1Space H1Space_from_c_H1Space(c_H1Space *h):
    cdef H1Space n
    n = <H1Space>PY_NEW(H1Space)
    n.thisptr = h
    return n

cdef Mesh Mesh_from_c_Mesh(c_Mesh *h):
    cdef Mesh n
    n = <Mesh>PY_NEW(Mesh)
    n.thisptr = h
    return n

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

    def h1_norm(self):
        cdef c_MeshFunction *m = <c_MeshFunction *>(self.thisptr)
        return h1_norm(m)

    def get_pt_value(self, x, y):
        cdef c_MeshFunction *m = <c_MeshFunction *>(self.thisptr)
        return m.get_pt_value(x, y)

    def get_mesh(self):
        cdef c_MeshFunction *m = <c_MeshFunction *>(self.thisptr)
        return Mesh_from_c_Mesh(m.get_mesh())

cdef class Solution(MeshFunction):

    def __cinit__(self):
        self.thisptr = <c_Function *>new_Solution()

    def __dealloc__(self):
        delete(self.thisptr)

    def set_zero(self, Mesh m):
        (<c_Solution *>(self.thisptr)).set_zero(m.thisptr)

    def set_fe_solution(self, H1Space s, PrecalcShapeset pss, ndarray v):
        """
        Sets the solution using the coefficient vector Y.
        """
        cdef int n = len(v)
        cdef int i
        from numpy import array
        cdef ndarray vec = array(v, dtype="double")
        cdef scalar *pvec = <scalar *>vec.data
        (<c_Solution *>(self.thisptr)).set_fe_solution(s.thisptr, pss.thisptr,
                pvec)

    # the get_fe_solution() method is is not yet implemented in the C++ hermes:
    #def get_fe_solution(self):
    #    """
    #    Returns the Y coefficients vector as a numpy array.
    #    """
    #    cdef int Ylen
    #    cdef scalar *Y
    #    (<c_Solution *>(self.thisptr)).get_fe_solution(&Ylen, &Y)
    #    if not (Ylen > 0 and Y != NULL):
    #        raise Exception("Ylen (%d) or Y is not valid." % Ylen)
    #    from numpy import empty
    #    cdef ndarray vec = empty([Ylen], dtype="double")
    #    cdef scalar *pvec = <scalar *>vec.data
    #    memcpy(pvec, Y, Ylen*sizeof(scalar))
    #    return vec

cdef class Filter(MeshFunction):
    pass

cdef class SimpleFilter(Filter):
    pass

cdef class VonMisesFilter(Filter):

    def __cinit__(self, MeshFunction sln1, MeshFunction sln2, double l, double m):
        self.thisptr = <c_Function *>new_VonMisesFilter(<c_MeshFunction *>sln1.thisptr, <c_MeshFunction *>sln2.thisptr, l, m)

    #def __dealloc__(self):
    #    delete(self.thisptr)

cdef class MagFilter(SimpleFilter):

    def __cinit__(self, MeshFunction sln1, MeshFunction sln2, int item1, int item2):
        self.thisptr = <c_Function *>new_MagFilter(<c_MeshFunction *>sln1.thisptr, <c_MeshFunction *>sln2.thisptr, item1, item2)

    #def __dealloc__(self):
    #    delete(self.thisptr)

cdef class DiffFilter(SimpleFilter):

    def __cinit__(self, MeshFunction sln1, MeshFunction sln2,
            int item1=FN_VAL, int item2=FN_VAL):
        self.thisptr = <c_Function *>new_DiffFilter(<c_MeshFunction *>sln1.thisptr, <c_MeshFunction *>sln2.thisptr, item1, item2)

    #def __dealloc__(self):
    #    delete(self.thisptr)

cdef class SumFilter(SimpleFilter):

    def __cinit__(self, MeshFunction sln1, MeshFunction sln2,
            int item1=FN_VAL, int item2=FN_VAL):
        self.thisptr = <c_Function *>new_SumFilter(<c_MeshFunction *>sln1.thisptr, <c_MeshFunction *>sln2.thisptr, item1, item2)

    #def __dealloc__(self):
    #    delete(self.thisptr)

cdef class SquareFilter(SimpleFilter):

    def __cinit__(self, MeshFunction sln, int item=FN_VAL):
        self.thisptr = <c_Function *>new_SquareFilter(<c_MeshFunction *>sln.thisptr, item)

    #def __dealloc__(self):
    #    delete(self.thisptr)

cdef class WeakForm:

    def __cinit__(self, int neq=1):
        self.thisptr = new_WeakForm(neq)

    def __dealloc__(self):
        delete(self.thisptr)

cdef class DummySolver(Solver):

    def __cinit__(self):
        self.thisptr = <c_Solver *>(new_DummySolver())

    def __dealloc__(self):
        delete(self.thisptr)

cdef class LinSystem:

    def __init__(self, WeakForm wf, Solver solver):
        self.thisptr = new_LinSystem(wf.thisptr, solver.thisptr)

    #def __dealloc__(self):
    #    delete(self.thisptr)

    def set_spaces(self, *args):
        self._spaces = args
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

    def set_pss(self, *args):
        self._pss = args
        cdef int n = len(args)
        cdef PrecalcShapeset s
        if n == 1:
            s = args[0]
            self.thisptr.set_pss(n, s.thisptr)
        else:
            raise NotImplementedError()

    def solve_system(self, *args):
        cdef int n = len(args)
        cdef Solution a, b, c
        cdef ndarray vec
        cdef scalar *pvec
        if n == 1:
            a = args[0]
            #self.thisptr.solve(n, a.thisptr)
            A = self.get_matrix()
            rhs = self.get_rhs()
            from scipy.sparse.linalg import cg
            x, res = cg(A, rhs)
            from numpy import array
            vec = array(x, dtype="double")
            pvec = <scalar *>vec.data
            (<c_Solution *>(a.thisptr)).set_fe_solution(
                    self.thisptr.get_space(0),
                    self.thisptr.get_pss(0),
                    pvec)
        elif n == 2:
            a, b = args
            self.thisptr.solve(n, a.thisptr, b.thisptr)
        elif n == 3:
            a, b, c = args
            self.thisptr.solve(n, a.thisptr, b.thisptr, c.thisptr)
        else:
            raise NotImplementedError()

    def assemble(self):
        self.thisptr.assemble()

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

    def get_rhs(self):
        """
        Return the RHS as a numpy array
        """
        cdef scalar *rhs
        cdef int n
        self.thisptr.get_rhs(rhs, n)
        from numpy import empty
        cdef ndarray vec = empty([n], dtype="double")
        cdef double *pvec = <double *>vec.data
        memcpy(pvec, rhs, n*sizeof(double))
        return vec

cdef class RefSystem(LinSystem):

    def __init__(self, LinSystem ls):
        self.thisptr = <c_LinSystem *>new_RefSystem(ls.thisptr)

    def assemble(self):
        (<c_RefSystem *>(self.thisptr)).assemble()

    # this is commented out, because get_ref_space() is not yet implemented in
    # C++ hermes2d
    #def get_ref_space(self, int eq):
    #    cdef c_H1Space *r = <c_H1Space *>(
    #            (<c_RefSystem *>(self.thisptr)).get_ref_space(eq)
    #        )
    #    return H1Space_from_c_H1Space(r)


#cdef class DiscreteProblem:
#
#    def __cinit__(self):
#        self.thisptr = new_DiscreteProblem()
#
#    def __dealloc__(self):
#        delete(self.thisptr)
#
#    def copy(self, DiscreteProblem ep):
#        self.thisptr.copy(ep.thisptr)
#
#    def set_num_equations(self, int neq):
#        self.thisptr.set_num_equations(neq)
#
#    def set_spaces(self, *args):
#        cdef int n = len(args)
#        cdef H1Space a, b, c
#        if n == 1:
#            a = args[0]
#            self.thisptr.set_spaces(n, a.thisptr)
#        elif n == 2:
#            a, b = args
#            self.thisptr.set_spaces(n, a.thisptr, b.thisptr)
#        elif n == 3:
#            a, b, c = args
#            self.thisptr.set_spaces(n, a.thisptr, b.thisptr, c.thisptr)
#        else:
#            raise NotImplementedError()
#
#    def set_external_fns(self, *args):
#        cdef int n = len(args)
#        cdef Solution a, b, c
#        if n == 1:
#            a = args[0]
#            self.thisptr.set_external_fns(n, a.thisptr)
#        elif n == 2:
#            a, b = args
#            self.thisptr.set_external_fns(n, a.thisptr, b.thisptr)
#        elif n == 3:
#            a, b, c = args
#            self.thisptr.set_external_fns(n, a.thisptr, b.thisptr, c.thisptr)
#        else:
#            raise NotImplementedError()
#
#    def set_pss(self, *args):
#        cdef int n = len(args)
#        cdef PrecalcShapeset s
#        if n == 1:
#            s = args[0]
#            self.thisptr.set_pss(n, s.thisptr)
#        else:
#            raise NotImplementedError()
#
#    #def set_bilinear_form(self, int i, int j, BiFormFnVol unsym):
#    #    self.thisptr.set_bilinear_form(i, j, unsym)
#
#    def create_matrix(self):
#        self.thisptr.create_matrix()
#
#    def set_quiet(self, quiet=True):
#        self.thisptr.set_quiet(quiet)
#
#    def assemble_matrix_and_rhs(self):
#        self.thisptr.assemble_matrix_and_rhs()
#
#    def solve_system(self, *args):
#        cdef int n = len(args)
#        cdef Solution a, b, c
#        if n == 1:
#            a = args[0]
#            self.thisptr.solve_system(n, a.thisptr)
#        elif n == 2:
#            a, b = args
#            self.thisptr.solve_system(n, a.thisptr, b.thisptr)
#        elif n == 3:
#            a, b, c = args
#            self.thisptr.solve_system(n, a.thisptr, b.thisptr, c.thisptr)
#        else:
#            raise NotImplementedError()
#
#    def save_matrix_matlab(self, char *filename, char *varname):
#        self.thisptr.save_matrix_matlab(filename, varname)
#
#    def save_matrix_coo(self, char *filename):
#        self.thisptr.save_matrix_coo(filename)
#
#    def get_matrix_csc(self):
#        """
#        Returns the matrix A as a (Ap, Ai, Ax) tuple.
#
#        See also DiscreteProblem.get_matrix() to get a SciPy matrix.
#        """
#        cdef int *Ap, *Ai, n, nnz
#        cdef scalar *Ax
#        self.thisptr.get_matrix(Ap, Ai, Ax, n)
#        nnz = Ap[n]
#        from numpy import zeros
#        cdef ndarray aAp = zeros(n+1, dtype="int32")
#        cdef ndarray aAi = zeros(nnz, dtype="int32")
#        cdef ndarray aAx = zeros(nnz, dtype="double")
#        cdef int *pAp = <int *>aAp.data
#        cdef int *pAi = <int *>aAi.data
#        cdef double *pAx = <double *>aAx.data
#        cdef int i
#        for i in range(n+1):
#            pAp[i] = Ap[i]
#        for i in range(nnz):
#            pAi[i] = Ai[i]
#            pAx[i] = Ax[i]
#        return aAp, aAi, aAx
#
#    def get_matrix(self):
#        """
#        Returns the global matrix A as a SciPy matrix.
#        """
#        from scipy.sparse import csc_matrix
#        Ap, Ai, Ax = self.get_matrix_csc()
#        return csc_matrix((Ax, Ai, Ap))
#
#    def free_matrix(self):
#        self.thisptr.free_matrix()

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

    def __dealloc__(self):
        delete(self.thisptr)

    def calc_error(self, MeshFunction sln, MeshFunction rsln):
        return self.thisptr.calc_error(sln.thisptr, rsln.thisptr)

    def adapt(self, double thr, int strat=0, h_only=False):
        self.thisptr.adapt(thr, strat, h_only)

cdef class H1OrthoHP:
    cdef c_H1OrthoHP *thisptr

    def __cinit__(self, *args):
        cdef int n = len(args)
        cdef H1Space a, b, c, d
        if n == 1:
            a = args[0]
            self.thisptr = new_H1OrthoHP(n, a.thisptr)
        elif n == 2:
            a, b = args
            self.thisptr = new_H1OrthoHP(n, a.thisptr, b.thisptr)
        elif n == 3:
            a, b, c = args
            self.thisptr = new_H1OrthoHP(n, a.thisptr, b.thisptr, c.thisptr)
        elif n == 4:
            a, b, c, d = args
            self.thisptr = new_H1OrthoHP(n, a.thisptr, b.thisptr, c.thisptr,
                    d.thisptr)
        else:
            raise NotImplementedError()

    def __dealloc__(self):
        delete(self.thisptr)

    def calc_error(self, MeshFunction sln, MeshFunction rsln):
        return self.thisptr.calc_error(sln.thisptr, rsln.thisptr)

    def calc_error_4(self, sln_list, rsln_list):
        return self.thisptr.calc_error_n(4,
                (<MeshFunction>(sln_list[0])).thisptr,
                (<MeshFunction>(sln_list[1])).thisptr,
                (<MeshFunction>(sln_list[2])).thisptr,
                (<MeshFunction>(sln_list[3])).thisptr,
                (<MeshFunction>(rsln_list[0])).thisptr,
                (<MeshFunction>(rsln_list[1])).thisptr,
                (<MeshFunction>(rsln_list[2])).thisptr,
                (<MeshFunction>(rsln_list[3])).thisptr,
                )

    def adapt(self, double thr, int strat=0, h_only=False):
        self.thisptr.adapt(thr, strat, h_only)

cdef class Linearizer:
    """
    Linearizes the solution.

    It returns the triangles and vertices and you can then use it to visualize
    the solution.

    Example:

        In [40]: l = Linearizer(sln)

        In [44]: l.process_solution(sln)
        Linearizer: 5519 verts, 10847 tris in 0.05 sec

        In [45]: l.get_vertices()
        Out[45]:
        array([[  0.00000000e+00,  -1.00000000e+00,  -2.22396971e-17],
               [  1.00000000e+00,  -1.00000000e+00,  -1.64798730e-17],
               [ -1.00000000e+00,   0.00000000e+00,   8.09899023e-17],
               ...,
               [  1.48437500e-01,  -1.56250000e-02,   1.62359362e-01],
               [  1.32812500e-01,   0.00000000e+00,   1.56012622e-01],
               [  1.32812500e-01,  -1.56250000e-02,   1.50562411e-01]])

        In [46]: l.get_triangles()
        Out[46]:
        array([[   3, 5448,   29],
               [  27, 5445,   28],
               [  29,   28,   26],
               ...,
               [5499, 5498, 5479],
               [5510, 5493, 5491],
               [5513, 5508, 5491]], dtype=int32)

    """

    cdef c_Linearizer *thisptr

    def __cinit__(self):
        self.thisptr = new_Linearizer()

    def __dealloc__(self):
        delete(self.thisptr)

    def process_solution(self, MeshFunction sln):
        cdef Function s = <Function>sln
        self.thisptr.process_solution(<c_MeshFunction *>sln.thisptr)

    def get_vertices(self):
        """
        Returns the list of vertices.

        It's a list of triples, where each triple contains (x, y, val), where
        x, y are the "x" and "y" coordinates of the vertex and "val" is the
        value of the solution at the vertex.

        Example:

        In [45]: l.get_vertices()
        Out[45]:
        array([[  0.00000000e+00,  -1.00000000e+00,  -2.22396971e-17],
               [  1.00000000e+00,  -1.00000000e+00,  -1.64798730e-17],
               [ -1.00000000e+00,   0.00000000e+00,   8.09899023e-17],
               ...,
               [  1.48437500e-01,  -1.56250000e-02,   1.62359362e-01],
               [  1.32812500e-01,   0.00000000e+00,   1.56012622e-01],
               [  1.32812500e-01,  -1.56250000e-02,   1.50562411e-01]])
        """
        cdef double3 *vert = self.thisptr.get_vertices()
        cdef int nvert = self.thisptr.get_num_vertices()
        from numpy import empty
        cdef ndarray vec = empty([3*nvert], dtype="double")
        cdef double *pvec = <double *>vec.data
        memcpy(pvec, vert, 3*nvert*sizeof(double))
        return vec.reshape((nvert, 3))

    def get_num_vertices(self):
        return self.thisptr.get_num_vertices()

    def get_triangles(self):
        """
        Returns a list of triangles.

        The list contains triples of vertices IDs. Use get_vertices() to obtain
        vertices coordinates.

        Example:

        In [46]: l.get_triangles()
        Out[46]:
        array([[   3, 5448,   29],
               [  27, 5445,   28],
               [  29,   28,   26],
               ...,
               [5499, 5498, 5479],
               [5510, 5493, 5491],
               [5513, 5508, 5491]], dtype=int32)
        """
        cdef int3 *tri = self.thisptr.get_triangles()
        cdef int ntri = self.thisptr.get_num_triangles()
        from numpy import empty
        cdef ndarray vec = empty([3*ntri], dtype="int32")
        cdef int *pvec = <int *>vec.data
        memcpy(pvec, tri, 3*ntri*sizeof(int))
        return vec.reshape((ntri, 3))

    def get_num_triangles(self):
        return self.thisptr.get_num_triangles()

    def get_edges(self):
        """
        Returns a list of edges.

        The list contains triples of vertices IDs. Use get_vertices() to obtain
        vertices coordinates.

        Example:

        In [47]: l.get_edges()
        Out[47]:
        array([[   3,   27,    0],
               [  27,   24,    0],
               [  24,   30,    0],
               ...,
               [5339, 5070,    4],
               [5070, 5077,    4],
               [5077,   11,    4]], dtype=int32)
        """
        cdef int3 *edges = self.thisptr.get_edges()
        cdef int nedges = self.thisptr.get_num_edges()
        from numpy import empty
        cdef ndarray vec = empty([3*nedges], dtype="int32")
        cdef int *pvec = <int *>vec.data
        memcpy(pvec, edges, 3*nedges*sizeof(int))
        return vec.reshape((nedges, 3))

    def get_num_edges(self):
        return self.thisptr.get_num_edges()

    def get_min_value(self):
        return self.thisptr.get_min_value()

    def get_max_value(self):
        return self.thisptr.get_max_value()

    def save_data(self, char* filename):
        self.thisptr.save_data(filename)

    def load_data(self, char* filename):
        self.thisptr.load_data(filename)

cdef class View:

    def wait(self):
        View_wait()

cdef class ScalarView(View):
    cdef c_ScalarView *thisptr

    def __cinit__(self, char *title, int x=-1, int y=-1, int width=1000,
                int height = 800):
        self.thisptr = new_ScalarView(title, x, y, width, height)

    #def __dealloc__(self):
    #    delete(self.thisptr)

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

    def wait(self):
        self.thisptr.wait()

cdef class BaseView(View):
    cdef c_BaseView *thisptr

    def __cinit__(self, char *title="BaseView", int x=-1, int y=-1,
            int width=1000, int height = 800):
        self.thisptr = new_BaseView(title, x, y, width, height)

    #def __dealloc__(self):
    #    delete(self.thisptr)

    def show(self, H1Space s):
        self.thisptr.show(s.thisptr)

cdef class MeshView(View):
    cdef c_MeshView *thisptr

    def __cinit__(self, char *title="MeshView", int x=-1, int y=-1,
            int width=1000, int height = 800):
        self.thisptr = new_MeshView(title, x, y, width, height)

    #def __dealloc__(self):
    #    delete(self.thisptr)

    def show(self, Mesh s):
        self.thisptr.show(s.thisptr)

cdef class OrderView(View):
    cdef c_OrderView *thisptr

    def __cinit__(self, char *title="OrderView", int x=-1, int y=-1,
            int width=1000, int height = 800):
        self.thisptr = new_OrderView(title, x, y, width, height)

    #def __dealloc__(self):
    #    delete(self.thisptr)

    def show(self, H1Space s):
        self.thisptr.show(s.thisptr)

cdef class VectorView(View):
    cdef c_VectorView *thisptr

    def __cinit__(self, char *title, int x=-1, int y=-1, int width=1000,
            int height = 800):
        self.thisptr = new_VectorView(title, x, y, width, height)

    #def __dealloc__(self):
    #    delete(self.thisptr)

    def show(self, MeshFunction s1, MeshFunction s2, eps = 0.008):
        (<c_VectorView *>(self.thisptr)).show(<c_MeshFunction *>s1.thisptr, <c_MeshFunction *>s2.thisptr, eps)

    def show_scale(self, show=True):
        self.thisptr.show_scale(show)

    def set_min_max_range(self, double min, double max):
        self.thisptr.set_min_max_range(min, max)

#cdef class MatrixView(View):
#    cdef c_MatrixView *thisptr
#
#    def __cinit__(self, char *title, int x=-1, int y=-1, int width=1000,
#            int height = 800):
#        self.thisptr = new_MatrixView(title, x, y, width, height)

    #def __dealloc__(self):
    #    delete(self.thisptr)

    #def show(self, DiscreteProblem ep):
    #    self.thisptr.show(ep.thisptr)

def init_hermes2d_wrappers():
    init_global_empty_tuple()

def set_verbose(verbose):
    """
    Sets the global verbose_mode variable.

    That variable controls how verbose hermes2d is.

    Returns the old status.
    """
    global c_verbose_mode
    global c_info_mode
    old_flag = c_verbose_mode
    c_verbose_mode = verbose
    c_info_mode = verbose
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

#def glut_main_loop():
#    """
#    This function waits for all Views to finish.
#
#    E.g. until the user closes them. If you don't call this function, all
#    windows will be closed when the program ends, so you need to call it if you
#    want to give the user a chance to manipulate with the plots at the end of
#    your script (program).
#    """
#    finish_glut_main_loop(False)

init_hermes2d_wrappers()
