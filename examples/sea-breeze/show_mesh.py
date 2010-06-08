#! /usr/bin/env python

import sys
from hermes2d import Mesh

m = Mesh()
m.load(sys.argv[1])
m.plot(lib="mpl", method="orders")
