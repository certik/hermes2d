#! /usr/bin/env python

"""
This script derives the initial conditions for the sea breeze model.

"""
import os

from jinja2 import Environment, FileSystemLoader
from sympy import var, ccode

var("z")

# Everything is in SI units
g = 9.80665
R = 287.14
c_v = 20.8
p = (100000 - 11476.*(z) + 529.54*(z)*(z) - 9.38*(z)*(z)*(z)).subs(z, z/1000)
p_0 = p.subs(z, 0)

rho = -p.diff(z)/g
rho_0 = rho.subs(z, 0)
T = p/(rho*R)
T_0 = T.subs(z, 0)

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
