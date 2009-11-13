#! /usr/bin/env python

"""
This script derives the initial conditions for the sea breeze model.

"""
import os
from math import tanh

from jinja2 import Environment, FileSystemLoader
from sympy import var, ccode

var("z")

# Everything is in SI units

# physical constants:
g = 9.80665
R = 287.14
c_v = 20.8

# values of pressure and temperature at the sea level:
p_0 = 10**5
T_0 = 300.5

# O(z**2) terms in the pressure are given by empirical data:
# for now we prescribe 0, so that the equations are linear
#p_higher_terms = 0.00052954*z**2 - 9.38e-9*z**3
p_higher_terms = 0

# The rest is given by the euler equations and the ideal gas law:
rho_0 = p_0/(R*T_0)
p = (p_0 - rho_0*g*z + p_higher_terms)
rho = -p.diff(z)/g
T = p/(rho*R)
T = T.series(z, 0, 4).removeO()

print "Equations in SI units"
print
print "p =", p
print "rho =", rho
print "T =", T
print

def print_values_h(h):
    print "values at z = %dm:" % h
    print "p(%d) = %f" % (h, p.subs(z, h))
    print "rho(%d) = %f" % (h, rho.subs(z, h))
    print "T(%d) = %f" % (h, T.subs(z, h))
    print "E(%d) = %f" % (h, rho.subs(z, h) * T.subs(z, h) * c_v)

print_values_h(0)
print
print_values_h(4000)

if 0:
    print "-"*80
    print "Equations in SI units, 'z' in km"
    print
    print "p =", p.subs(z, z*1000)
    print "rho =", rho.subs(z, z*1000)
    print "T =", T.subs(z, z*1000)

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
