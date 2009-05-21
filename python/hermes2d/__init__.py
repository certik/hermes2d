# Just import everyting from our cython module:

import os

from _hermes2d import *
from plot import sln2png, plot_sln_mpl, plot_sln_mayavi, ScalarView
from runtests import test

def get_pxd_include():
    """
    Returns an absolute path to *.pxd files that are needed in order to build
    something against hermes2d.
    """
    return get_include()

def get_include():
    """
    Return the directory in the package that contains the hermes2d/*.h header
    files.

    Extension modules that need to compile against hermes2d should use this
    function to locate the appropriate include directory. Using distutils:

      import hermes2d
      Extension('extension_name', ...
                include_dirs=[hermes2d.get_include()])
    """
    this_dir = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(this_dir, "include")

def get_lib():
    """
    Returns the path to the *.so libraries that one needs to link
    against.
    """
    this_dir = os.path.abspath(os.path.dirname(__file__))
    return this_dir
