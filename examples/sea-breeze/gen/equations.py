#! /usr/bin/env python

"""
This script derives the initial conditions for the sea breeze model.

"""
import os

from jinja2 import Environment, FileSystemLoader
from sympy import var, ccode

var("z")

# Everything is in SI units

# physical constants:
g = 9.80665
R = 287.14
c_v = 20.8

# values of the various quantities at the sea-level:
p_0 = 10**5
T_0 = 297.602407347392

# One also prescribes the higher terms in the expansion of "p". The rest
# is given by the euler equations and the ideal gas law:
rho_0 = p_0/(R*T_0)
p = (p_0 - rho_0*g*z + 0.00052954*z**2 - 9.38e-9*z**3)
rho = -p.diff(z)/g
T = p/(rho*R)
T = T.series(z, 0, 4).removeO()

print "-"*80
print "Equations in SI units"
print
print "p =", p
print "rho =", rho
print "T =", T

def format_code(e):
    return "(" + ccode(e).replace("z", "(z)") + ")"

params = {
        "p_z": format_code(p),
        "rho_z": format_code(rho),
        "T_z": format_code(T),
        "R": R,
        "g": g,
        "c_v": c_v,
        }

template = "params.h"
env = Environment(loader=FileSystemLoader('.'))
t = env.get_template(template)
open(os.path.join("..", template), "w").write(t.render({
        "params": params,
        }))
