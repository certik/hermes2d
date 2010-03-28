from numpy import array
from numpy.linalg import norm

import pylab

from hermes2d import Mesh, set_verbose, MeshView

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

    def __plot__(self):
        styles = {
                0: {"color": "black", "lw": 1},
                1: {"color": "blue", "lw": 2},
                2: {"color": "green", "lw": 2},
                3: {"color": "red", "lw": 2},
                4: {"color": "orange", "lw": 2},
            }
        p0 = self.get_point_0()
        p1 = self.get_point_1()
        pylab.plot([p0[0], p1[0]], [p0[1], p1[1]], **styles[self.marker])
        # normal
        p0 = self.get_point_middle()
        d = self.normal * 0.1
        p1 = p1 + d
        #pylab.arrow(p0[0], p0[1], d[0], d[1], **styles[self.marker])
        pylab.gca().add_patch(pylab.Arrow(p0[0], p0[1], d[0], d[1],
            width = 0.05, **styles[self.marker]))

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

    def __plot__(self):
        for e in self._edges:
            self._edges[e].__plot__()
        pylab.gca().set_aspect("equal")

    def plot(self):
        self.__plot__()
        pylab.show()

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
    print edges
    print "Done."

    edges.plot()
    #mview = MeshView()
    #mview.show(mesh, lib="mpl", method="orders")

main()
