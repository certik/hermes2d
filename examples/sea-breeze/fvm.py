from math import atan2, sin, cos

from numpy import array, dot, zeros
from numpy.linalg import norm
import pylab

from hermes2d import Mesh, set_verbose, MeshView
from _numerical_flux import numerical_flux
from _forms import R as R_const, c_v

marker_bottom = 1
marker_right  = 2
marker_top    = 3
marker_left   = 4


class Plot(object):

    def __init__(self):
        self._num_markers = 5
        self._styles = {
                0: {"color": "black", "lw": 1},
                1: {"color": "blue", "lw": 2},
                2: {"color": "green", "lw": 2},
                3: {"color": "red", "lw": 2},
                4: {"color": "orange", "lw": 2},
            }
        self.color_count = [0]*self._num_markers

    def add_edge(self, p0, p1, normal, marker):
        assert marker < self._num_markers
        self.color_count[marker] += 1
        if self.color_count[marker] == 1:
            pylab.plot([p0[0], p1[0]], [p0[1], p1[1]], label=str(marker),
                    **self._styles[marker])
        else:
            pylab.plot([p0[0], p1[0]], [p0[1], p1[1]], **self._styles[marker])
        # normal
        p0 = (p0+p1)/2
        d = normal * 0.1
        pylab.gca().add_patch(pylab.Arrow(p0[0], p0[1], d[0], d[1],
            width = 0.05, color=self._styles[marker]["color"]))

    def show(self):
        pylab.gca().set_aspect("equal")
        pylab.legend()
        pylab.show()

class Edge(object):

    def __init__(self, i, marker, pair_nodes):
        self._pair = pair_nodes
        self._normal = self.calculate_normal()
        self._elements = [i]
        self._boundary = True
        self._marker = marker

    def set_second_element(self, i):
        self._elements.append(i)
        self._boundary = False

    def calculate_normal(self):
        p0 = self.get_point_0()
        p1 = self.get_point_1()
        t = array([p1[0]-p0[0], p1[1]-p0[1]])
        t = t/norm(t)
        return array([t[1], -t[0]])

    @property
    def boundary(self):
        return self._boundary

    @property
    def marker(self):
        return self._marker

    @property
    def normal(self):
        return self._normal

    def get_point_0(self):
        return array(self._pair[0].coord)

    def get_point_1(self):
        return array(self._pair[1].coord)

    def get_point_middle(self):
        return (self.get_point_0() + self.get_point_1()) / 2

    @property
    def elements(self):
        return self._elements

    def __str__(self):
        s = "elements: %s, boundary: %s, marker: %s, n=%s" % (self.elements,
                self.boundary, self.marker, self.normal)
        return s

    def __plot__(self, plot):
        p0 = self.get_point_0()
        p1 = self.get_point_1()
        plot.add_edge(p0, p1, self.normal, self.marker)

class Edges(object):

    def __init__(self, mesh):
        self._nodes_dict = mesh.nodes_dict
        self._elements = mesh.elements
        self._nodes = self.extract_nodes(self._elements)
        self._edges = self.extract_edges(mesh)

    def extract_nodes(self, elements):
        nodes = list(set(array(elements).flat))
        nodes.sort()

        _nodes = self._nodes_dict.keys()
        _nodes.sort()
        assert nodes == _nodes
        return nodes

    @property
    def edges(self):
        return self._edges

    def extract_edges(self, mesh):
        edges = {}
        for i in range(mesh.num_elements):
            e = mesh.get_element(i)
            if not e.active:
                continue
            nodes_edge = e.nodes_edge
            nodes_vertex = e.nodes_vertex
            for j, ed in enumerate(nodes_edge):
                jp1 = j + 1
                if jp1 >= len(nodes_edge):
                    jp1 = 0
                pair_nodes = (nodes_vertex[j], nodes_vertex[jp1])
                pair = (pair_nodes[0].id, pair_nodes[1].id)
                if pair in edges:
                    raise Exception("This edge was already processed")
                pair_reversed = (pair[1], pair[0])
                if pair_reversed in edges:
                    edges[pair_reversed].set_second_element(i)
                    continue
                edges[pair] = Edge(i, ed.marker, pair_nodes)
        return edges

    def __str__(self):
        s = "Edges:\n"
        for edge in self._edges:
            s += "%12s:  %s\n" % (edge, self._edges[edge])
        return s

    def __plot__(self, plot):
        for e in self._edges:
            self._edges[e].__plot__(plot)

    def plot(self):
        p = Plot()
        self.__plot__(p)
        p.show()

def T_rot(beta):
    # this is the 3D rotation matrix (2.2.51 in Pavel's master thesis)
    # in 2D without the middle column and row, alpha = 0
    alpha = 0
    return array([
        [1, 0, 0, 0],
        [0, cos(alpha)*cos(beta), sin(beta), 0],
        [0, -cos(alpha)*sin(beta), cos(beta), 0],
        [0, 0, 0, 1]
        ])

def calc_p(w):
    w0, w1, w3, w4 = w
    p = R_const/c_v * (w4 - (w1**2 + w3**2)/(2*w0))
    return p

def calculate_flux(edge, state_on_elements):
    w_l = state_on_elements[edge.elements[0]]
    if edge.boundary:
        if edge.marker in [marker_left, marker_right]:
            w_r = array([1., 50., 0., 1.e5])
        elif edge.marker in [marker_top, marker_bottom]:
            alpha = atan2(edge.normal[1], edge.normal[0])
            p = calc_p(w_l)
            flux_local = array([0., p, 0., 0.])
            return dot(T_rot(alpha), flux_local)
        else:
            raise Exception("Unhandled boundary")
    else:
        w_r = state_on_elements[edge.elements[1]]
    return numerical_flux(w_l, w_r, edge.normal)

def assembly(edges, state_on_elements):
    elem_contrib = {}
    for e in state_on_elements:
        elem_contrib[e] = zeros((4,))
    for e in edges.edges:
        edge = edges.edges[e]
        flux = calculate_flux(edge, state_on_elements)
        if edge.boundary:
            elem_contrib[edge.elements[0]] += flux
        else:
            elem_contrib[edge.elements[0]] += flux
            elem_contrib[edge.elements[1]] -= flux
    for e in elem_contrib:
        print "-"*80
        print e
        print elem_contrib[e]

def main():
    set_verbose(False)
    mesh = Mesh()
    print "Loading mesh..."
    mesh.load("GAMM-channel.mesh")
    mesh.refine_element(1, 2)
    #mesh.refine_all_elements()
    #mesh.refine_all_elements()
    #mesh.refine_all_elements()
    #mesh.refine_all_elements()

    print "Constructing edges..."
    nodes = mesh.nodes_dict
    edges = Edges(mesh)
    elements = mesh.elements
    print "Done."

    print "Assembly..."
    state_on_elements = {}
    for e in mesh.active_elements:
        state_on_elements[e.id] = array([1., 50., 0., 1.e5])
    assembly(edges, state_on_elements)
    print "Done."

    edges.plot()
    #mview = MeshView()
    #mview.show(mesh, lib="mpl", method="orders")

main()
