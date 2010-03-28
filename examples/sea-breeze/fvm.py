from hermes2d import Mesh, set_verbose, MeshView

set_verbose(False)
mesh = Mesh()
mesh.load("GAMM-channel.mesh")
mesh.refine_element(1, 2)
mesh.refine_all_elements()
mesh.refine_all_elements()
mesh.refine_all_elements()
mesh.refine_all_elements()

mview = MeshView()
mview.show(mesh, lib="mpl", method="orders")
