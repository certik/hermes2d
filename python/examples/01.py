#! /usr/bin/env python

from hermes2d import Mesh, MeshView

mesh = Mesh()
mesh.load("domain.mesh")
#mesh.refine_element(2)
mesh.refine_all_elements()
mesh.refine_all_elements()
mesh.refine_all_elements()

mview = MeshView("Hello world!", 100, 100, 500, 500)
mview.show(mesh)
