#! /usr/bin/env python

"""
This script derives the initial conditions for the sea breeze model.

"""
from sympy import var

var("z")

# Everything is in SI units
g = 9.80665
R = 287.14
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

print "-"*80
print "Equations in SI units, except 'z' in km"
print
print "p =", p.subs(z, z*1000)
print "rho =", rho.subs(z, z*1000)
print "T =", T.subs(z, z*1000)

print "-"*80
print "Values at z = 0m:"
print "p =", p.subs(z, 0)
print "rho =", rho.subs(z, 0)
print "T =", T.subs(z, 0)
print "Values at z = 11000m:"
print "p =", p.subs(z, 11000)
print "rho =", rho.subs(z, 11000)
print "T =", T.subs(z, 11000)
