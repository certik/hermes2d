import os

def get_example_mesh():
    """
    Returns an example mesh, that is distributed with hermes2d.
    """
    this_dir = os.path.dirname(os.path.abspath(__file__))
    mesh = os.path.join(this_dir, "domain.mesh")
    return os.path.normpath(mesh)

def get_sample_mesh():
    """
    Returns an example mesh, that is distributed with hermes2d.
    """
    this_dir = os.path.dirname(os.path.abspath(__file__))
    mesh = os.path.join(this_dir, "sample.mesh")
    return os.path.normpath(mesh)

def get_cylinder_mesh():
    """
    Returns an example mesh, that is distributed with hermes2d.
    """
    this_dir = os.path.dirname(os.path.abspath(__file__))
    mesh = os.path.join(this_dir, "cylinder4.mesh")
    return os.path.normpath(mesh)
