from numpy import array

from hermes2d import Mesh, set_verbose, MeshView

class Edge(object):

    def __init__(self, i, marker, pair):
        self._pair = pair
        self._elements = [i]
        self._boundary = True
        self._marker = marker

    def set_second_element(self, i):
        self._elements.append(i)
        self._boundary = False

    @property
    def boundary(self):
        return self._boundary

    @property
    def marker(self):
        return self._marker

    @property
    def elements(self):
        return self._elements

    def __str__(self):
        s = "elements: %s, boundary: %s, marker: %s" % (self.elements,
                self.boundary, self.marker)
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
                pair = (nodes_vertex[j].id, nodes_vertex[jp1].id)
                if pair in edges:
                    raise Exception("This edge was already processed")
                pair_reversed = (pair[1], pair[0])
                if pair_reversed in edges:
                    edges[pair_reversed].set_second_element(i)
                    continue
                edges[pair] = Edge(i, ed.marker, pair)
        return edges

    def __str__(self):
        s = "Edges:\n"
        for edge in self._edges:
            s += "%12s:  %s\n" % (edge, self._edges[edge])
        return s

def main():
    set_verbose(False)
    mesh = Mesh()
    mesh.load("GAMM-channel.mesh")
    mesh.refine_element(1, 2)

    nodes = mesh.nodes_dict
    edges = Edges(mesh)
    elements = mesh.elements
    print edges
    #mesh.refine_all_elements()
    #mesh.refine_all_elements()
    #mesh.refine_all_elements()
    #mesh.refine_all_elements()

    mview = MeshView()
    #mview.show(mesh, lib="mpl", method="orders")

main()
