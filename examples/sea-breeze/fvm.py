from numpy import array
from numpy.linalg import norm

from hermes2d import Mesh, set_verbose, MeshView

class Edge(object):

    def __init__(self, i, marker, pair_nodes):
        self._pair = pair_nodes
        self._normal = self.calculate_normal(pair_nodes)
        self._elements = [i]
        self._boundary = True
        self._marker = marker

    def set_second_element(self, i):
        self._elements.append(i)
        self._boundary = False

    def calculate_normal(self, pair_nodes):
        p0 = pair_nodes[0].coord
        p1 = pair_nodes[1].coord
        t = array([p1[0]-p0[0], p1[1]-p0[1]])
        t = t/norm(t)
        return (t[1], -t[0])

    @property
    def boundary(self):
        return self._boundary

    @property
    def marker(self):
        return self._marker

    @property
    def normal(self):
        return self._normal

    @property
    def elements(self):
        return self._elements

    def __str__(self):
        s = "elements: %s, boundary: %s, marker: %s, n=%s" % (self.elements,
                self.boundary, self.marker, self.normal)
        return s

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

    mview = MeshView()
    mview.show(mesh, lib="mpl", method="orders")

main()
