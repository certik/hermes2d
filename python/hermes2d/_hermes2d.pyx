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

cdef class Nurbs:
    cdef c_Nurbs *thisptr

    @property
    def degree(self):
        return self.thisptr.degree

    @property
    def np(self):
        return self.thisptr.np

    @property
    def pt(self):
        cdef double3 *pt = self.thisptr.pt
        cdef int np = self.thisptr.np
        cdef ndarray vec = array_double_c2numpy(<double *>pt, 3*np)
        return vec.reshape((np, 3))

    @property
    def nk(self):
        return self.thisptr.nk

    @property
    def kv(self):
        cdef double *kv = self.thisptr.kv
        cdef int nk = self.thisptr.nk
        cdef ndarray vec = array_double_c2numpy(<double *>kv, nk)
        return vec
        #return vec.reshape((kv,))

    @property
    def ref(self):
        return self.thisptr.nk

cdef class CurvMap:
    cdef c_CurvMap *thisptr

    @property
    def order(self):
        return self.thisptr.order

    @property
    def toplevel(self):
        return self.thisptr.toplevel

    def get_nurbs(self, int k):
        if self.thisptr.toplevel == 0:
            return None
        if self.thisptr.nurbs[k] == NULL:
            return None

        cdef Nurbs n = Nurbs()
        n.thisptr = self.thisptr.nurbs[k]
        return n

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
    def nodes_vertex_id(self):
        return [node.id for node in self.nodes_vertex]

    @property
    def nodes_edge(self):
        return [self.get_edge_node(i) for i in range(self.nvert)]

    def get_vertex_node(self, int id):
        cdef Node n = Node()
        n.thisptr = self.thisptr.vn[id]
        return n

    @property
    def curved_map(self):
        if self.thisptr.cm == NULL:
            return None

        cdef CurvMap cm = CurvMap()
        cm.thisptr = self.thisptr.cm
        return cm

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

def get_node_id(node):
    # This function is used in nodes.sort() in the Mesh class.
    # Cython doesn't yet support lambda functions, nor closures.
    return node.id


cdef class Mesh:
    cdef c_Mesh *thisptr

    def __cinit__(self):
        self.thisptr = new_Mesh()

    def __dealloc__(self):
        delete(self.thisptr)

    def copy(self, Mesh m):
        self.thisptr.copy(m.thisptr)

    def load(self, char* filename):
        cdef c_H2DReader mloader
        mloader.load(filename, self.thisptr)

    def load_str(self, char* mesh):
        cdef c_H2DReader mloader
        mloader.load_str(mesh, self.thisptr)

    @property
    def elements_markers(self):
        """
        Return the list of elements (as ids), together with their markers.

        This is equaivalent to the format of the elements argument to the
        Mesh.create() method.
        """
        element_list = []
        for i in range(self.num_elements):
            el = self.get_element(i)
            element_list.append(el.nodes_vertex_id + [el.marker])
        return element_list

    @property
    def elements(self):
        """
        Return the list of elements (as ids).
        """
        element_list = []
        for i in range(self.num_elements):
            el = self.get_element(i)
            if el.active:
                element_list.append(el.nodes_vertex_id)
        return element_list

    @property
    def nodes(self):
        """
        Returns a list of nodes coordinates.

        However, this is only usable if you don't refine.

        Use nodes_dict for more usable implementation.
        """
        # This is a really slow implementation, but it gets the job
        # done for now. Later on, we should get the list of nodes from C++
        # directly.
        nodes = []
        for i in range(self.num_elements):
            el = self.get_element(i)
            nodes.extend(el.nodes_vertex)
        node_dict = {}
        for node in nodes:
            node_dict[node.id] = node
        nodes = node_dict.values()
        nodes.sort(key=get_node_id)
        nodes_coord = [node.coord for node in nodes]
        return nodes_coord

    @property
    def nodes_dict(self):
        """
        Returns a dict of all active nodes coordinates.
        """
        nodes = []
        for i in range(self.num_elements):
            el = self.get_element(i)
            if el.active:
                nodes.extend(el.nodes_vertex)
        node_dict = {}
        for node in nodes:
            node_dict[node.id] = node.coord
        return node_dict

    @property
    def curves(self):
        """
        Return curves
        """
        crv = {}
        for i in range(self.num_elements):
            el = self.get_element(i)
            if el.active:
                for j in range(4):
                    cm = el.curved_map
                    if cm != None:
                        n = cm.get_nurbs(j)
                        if n != None:
                            sp = (n.pt[0][0], n.pt[0][1])
                            cp = (n.pt[1][0], n.pt[1][1])
                            ep = (n.pt[2][0], n.pt[2][1])
                            crv[i] = []
                            crv[i].append([sp, cp, ep])
        return crv

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

    def load_hermes(self, char* filename):
        """
        Loads a mesh (in the hermes format) from a file.

        This uses a pure Python reader. If you want to use the C++ flex reader,
        use load().
        """
        from mesh import read_hermes_format
        nodes, elements, boundary, nurbs = \
                read_hermes_format(filename)
        self.create(nodes, elements, boundary, nurbs)

    def get_elements_order(self, space):
        """
        Returns list of orders
        """
        orders_list = []
        for i in range(self.num_elements):
            el = self.get_element(i)
            if el.active:
                order = space.get_element_order(i)
                h = order & ((1 << 5) - 1)
                v = order >> 5

                import math
                ord = int(((h+v)/2.0))
                if ord == 0:
                    ord = 1

                #orders_list.append(int(((h+v)/2.0)))
                orders_list.append(ord)

        return orders_list

    def create(self, nodes, elements, boundary, nurbs):
        """
        Creates a mesh from a list of nodes, elements, boundary and nurbs.

        Example:

        >>> import hermes2d
        >>> m = hermes2d.Mesh()
        >>> m.create([
        ...         [0, -1],
        ...         [1, -1],
        ...         [-1, 0],
        ...         [0, 0],
        ...         [1, 0],
        ...         [-1, 1],
        ...         [0, 1],
        ...         [0.707106781, 0.707106781],
        ...     ], [
        ...         [0, 1, 4, 3, 0],
        ...         [3, 4, 7, 0],
        ...         [3, 7, 6, 0],
        ...         [2, 3, 6, 5, 0],
        ...     ], [
        ...         [0, 1, 1],
        ...         [1, 4, 2],
        ...         [3, 0, 4],
        ...         [4, 7, 2],
        ...         [7, 6, 2],
        ...         [2, 3, 4],
        ...         [6, 5, 2],
        ...         [5, 2, 3],
        ...     ], [
        ...         [4, 7, 45],
        ...         [7, 6, 45],
        ...     ])

        """
        if boundary is None:
            boundary = []
        if nurbs is None:
            nurbs = []
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

    def refine_towards_vertex(self, int marker, int depth):
        self.thisptr.refine_towards_vertex(marker, depth)

    def get_element(self, int id):
        cdef Element e = Element()
        e.thisptr = self.thisptr.get_element(id)
        return e

    def plot(self, *args, **kwargs):
        """
        Plots the mesh and shows it to the user.

        It passes all arguments to the MeshView.show() function, so read its
        documentation for the meaning.
        """
        from hermes2d import MeshView
        m = MeshView()
        m.show(self, *args, **kwargs)

    def convert_triangles_to_quads(self):
        self.thisptr.convert_triangles_to_quads()

    def convert_quads_to_triangles(self):
        self.thisptr.convert_quads_to_triangles()

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

    def get_element_order(self, int el_id):
        return self.thisptr.get_element_order(el_id)

    def get_num_dofs(self):
        return self.thisptr.get_num_dofs()

cdef api object H1Space_from_C(c_H1Space *h):
    cdef H1Space n
    n = <H1Space>PY_NEW(H1Space)
    n.thisptr = h
    return n

cdef api object Mesh_from_C(c_Mesh *h):
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
        return Mesh_from_C(m.get_mesh())

cdef api object Solution_from_C(c_Solution *s):
    cdef Solution n
    n = <Solution>PY_NEW(Solution)
    n.thisptr = <c_Function *>s
    return n

cdef class Solution(MeshFunction):

    def __cinit__(self):
        self.thisptr = <c_Function *>new_Solution()

    #def __dealloc__(self):
    #    delete(self.thisptr)

    def set_zero(self, Mesh m):
        (<c_Solution *>(self.thisptr)).set_zero(m.thisptr)

    def set_const(self, Mesh m, scalar):
        (<c_Solution *>(self.thisptr)).set_const(m.thisptr, scalar)


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
        cdef PrecalcShapeset s1, s2
        if n == 1:
            s1 = args[0]
            self.thisptr.set_pss(n, s1.thisptr)
        elif n == 2:
            s1 = args[0]
            s2 = args[1]
            self.thisptr.set_pss(n, s1.thisptr, s2.thisptr)
        else:
            raise NotImplementedError()

    def solve_system(self, *args, lib="scipy"):
        """
        Solves the linear system using scipy.
        """
        cdef int n = len(args)

        cdef Solution s0, s1, s2, s3
        cdef ndarray vec
        cdef scalar *pvec

        if lib == "hermes":
            if n == 1:
                s0 = args[0]
                self.thisptr.solve(n, s0.thisptr)
            elif n == 2:
                s0, s1 = args
                self.thisptr.solve(n, s0.thisptr, s1.thisptr)
            elif n == 3:
                s0, s1, s2 = args
                self.thisptr.solve(n, s0.thisptr, s1.thisptr, s2.thisptr)
            elif n == 4:
                s0, s1, s2, s3 = args
                self.thisptr.solve(n, s0.thisptr, s1.thisptr, s2.thisptr,
                        s3.thisptr)
            else:
                raise NotImplementedError()

        elif lib == "scipy":
            from scipy.sparse.linalg import cg
            from scipy.sparse.linalg import spsolve
            from numpy import array
            A = self.get_matrix()
            rhs = self.get_rhs()
            #x, res = cg(A, rhs)
            x = spsolve(A, rhs)
            vec = array(x, dtype="double")
            pvec = <scalar *>vec.data

            for i, sln in enumerate(args):
                (<c_Solution *>((<Solution>sln).thisptr)).set_fe_solution(
                        self.thisptr.get_space(i),
                        self.thisptr.get_pss(i),
                        pvec)
        else:
            raise NotImplementedError("Unknown library")

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
        aAp = array_int_c2numpy(Ap, n+1)
        aAi = array_int_c2numpy(Ai, nnz)
        aAx = array_double_c2numpy(Ax, nnz)
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
        return array_double_c2numpy(rhs, n)

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
    #    return H1Space_from_C(r)


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
    #cdef c_H1OrthoHP *thisptr

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
            self.thisptr = new_H1OrthoHP(n, a.thisptr, b.thisptr, c.thisptr, d.thisptr)
        else:
            raise NotImplementedError()

    def __dealloc__(self):
        delete(self.thisptr)


    def calc_error(self, MeshFunction sln, MeshFunction rsln):
        return self.thisptr.calc_error(sln.thisptr, rsln.thisptr)

    def calc_error_2(self, MeshFunction sln1, MeshFunction sln2, MeshFunction rsln1, MeshFunction rsln2):
        return self.thisptr.calc_error_2(sln1.thisptr, sln2.thisptr, rsln1.thisptr, rsln2.thisptr)

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

    Example::

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

    def __cinit__(self):
        self.thisptr = new_Linearizer()

    def __dealloc__(self):
        delete(self.thisptr)

    def process_solution(self, MeshFunction sln):
        self.thisptr.process_solution(<c_MeshFunction *>sln.thisptr)

    def get_vertices(self):
        """
        Returns the list of vertices.

        It's a list of triples, where each triple contains (x, y, val), where
        x, y are the "x" and "y" coordinates of the vertex and "val" is the
        value of the solution at the vertex.

        Example::

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
        cdef ndarray vec = array_double_c2numpy(<double *>vert, 3*nvert)
        return vec.reshape((nvert, 3))

    def get_num_vertices(self):
        return self.thisptr.get_num_vertices()

    def get_triangles(self):
        """
        Returns a list of triangles.

        The list contains triples of vertices IDs. Use get_vertices() to obtain
        vertices coordinates.

        Example::

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
        cdef ndarray vec = array_int_c2numpy(<int *>tri, 3*ntri)
        return vec.reshape((ntri, 3))

    def get_num_triangles(self):
        return self.thisptr.get_num_triangles()

    def get_edges(self):
        """
        Returns a list of edges.

        The list contains triples of vertices IDs. Use get_vertices() to obtain
        vertices coordinates.

        Example::

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
        cdef ndarray vec = array_int_c2numpy(<int *>edges, 3*nedges)
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

cdef class Vectorizer(Linearizer):

    def __cinit__(self):
        self.thisptr = <c_Linearizer *>new_Vectorizer()

    # this segfaults:
    #def __dealloc__(self):
    #    delete(<c_Vectorizer *>(self.thisptr))

    def process_solution(self, MeshFunction xsln, MeshFunction ysln,
            int xitem=c_FN_VAL_0, int yitem=c_FN_VAL_0, double eps=c_EPS_LOW):
        (<c_Vectorizer *>(self.thisptr)).process_solution(
                <c_MeshFunction *>xsln.thisptr, xitem,
                <c_MeshFunction *>ysln.thisptr, yitem,
                eps)

    def get_vertices(self):
        cdef c_Vectorizer *_self = <c_Vectorizer *>(self.thisptr)
        cdef double4 *vert = _self.get_vertices()
        cdef int nvert = self.thisptr.get_num_vertices()
        cdef ndarray vec = array_double_c2numpy(<double *>vert, 4*nvert)
        return vec.reshape((nvert, 4))

    def get_dashes(self):
        cdef c_Vectorizer *_self = <c_Vectorizer *>(self.thisptr)
        cdef int2 *dashes = _self.get_dashes()
        cdef int ndashes = _self.get_num_dashes()
        cdef ndarray vec = array_int_c2numpy(<int *>dashes, 2*ndashes)
        return vec.reshape((ndashes, 2))

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

import sys
import traceback

global_namespace = {"verbose": False}

cdef api void cmd(char *text):
    n = run_cmd(text, global_namespace)
    global_namespace.update(n)

cdef api void set_verbose_cmd(int verbose):
    global_namespace["verbose"] = verbose

cdef api void insert_int(char *name, int i):
    """
    Inserts the int "i" into the global namespace.

    Example:

    insert_int("a", 34);
    cmd("print a");

    This prints "34".
    """
    global_namespace.update({name: i})

cdef api void insert_double_array(char *name, double *A, int len):
    """
    Inserts an array of doubles into the global namespace as a NumPy array.

    Example:

    double a[3] = {1, 5, 3};
    insert_double_array("A", a, 3);
    cmd("print A");

    This prints "[ 1.  5.  3.]" (this is how the NumPy array is printed).
    """
    global_namespace.update({name: array_double_c2numpy(A, len)})

cdef api void insert_int_array(char *name, int *A, int len):
    """
    Inserts an array of ints into the global namespace as a NumPy array.

    Example:

    int a[3] = {1, 5, 3};
    insert_double_array("A", a, 3);
    cmd("print A");

    This prints "[1  5  3]" (this is how the NumPy array is printed).
    """
    global_namespace.update({name: array_int_c2numpy(A, len)})

cdef api void insert_object(char *name, object o):
    global_namespace.update({name: o})

cdef api object get_symbol(char *name):
    return global_namespace.get(name)

cdef ndarray array_int_c2numpy(int *A, int len):
    from numpy import empty
    cdef ndarray vec = empty([len], dtype="int32")
    cdef int *pvec = <int *>vec.data
    memcpy(pvec, A, len*sizeof(int))
    return vec

cdef ndarray array_double_c2numpy(double *A, int len):
    from numpy import empty
    cdef ndarray vec = empty([len], dtype="double")
    cdef double *pvec = <double *>vec.data
    memcpy(pvec, A, len*sizeof(double))
    return vec

cdef api void array_double_numpy2c_inplace(object A_n, double **A_c, int *n):
    cdef ndarray A = A_n
    if not (A.nd == 1 and A.strides[0] == sizeof(double)):
        from numpy import array
        A = array(A.flat, dtype="double")
    n[0] = len(A)
    A_c[0] = <double *>(A.data)

cdef api object run_cmd(char *text, object namespace):
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
        return namespace
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
