#! /usr/bin/env python

from hermes2d import Mesh, H1Shapeset, PrecalcShapeset, H1Space
from _forms import register_bc

mesh = Mesh()
mesh.load("domain-quad.mesh")

mesh.refine_all_elements()
mesh.refine_all_elements()
mesh.refine_all_elements()

shapeset_h1 = H1Shapeset()
pss_h1 = PrecalcShapeset(shapeset_h1)
s0 = H1Space(mesh, shapeset_h1)
s1 = H1Space(mesh, shapeset_h1)
s3 = H1Space(mesh, shapeset_h1)
s4 = H1Space(mesh, shapeset_h1)

register_bc(s0, s1, s3, s4)
