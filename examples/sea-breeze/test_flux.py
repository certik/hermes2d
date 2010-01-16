from math import sin, cos, pi

from numpy import array, zeros, dot
from numpy.linalg import inv

from _numerical_flux import matrix_R, matrix_R_inv, matrix_D_minus, \
        c_v
from _numerical_flux import R as R_const

def R(w):
    A = zeros((4, 4))
    for i in range(4):
        for j in range(4):
            A[i, j] = matrix_R(i, j, *w)
    return A

def R_inv(w):
    A = zeros((4, 4))
    for i in range(4):
        for j in range(4):
            A[i, j] = matrix_R_inv(i, j, *w)
    return A

def D_minus(w):
    A = zeros((4, 4))
    for i in range(4):
        for j in range(4):
            A[i, j] = matrix_D_minus(i, j, *w)
    return A

def f_x(w):
    w0, w1, w3, w4 = w
    p = R_const/c_v * (w4 - (w1**2 + w3**2)/(2*w0))
    return array([w1, w1**2/w0 + p, w1*w3/w0, w1/w0 * (w4 + p)])

def A_minus(w):
    return dot(R(w), dot(D_minus(w), R_inv(w)))

def f_riemann(w_l, w_r):
    return f_x(w_l) + dot(A_minus(w_r), w_r) - dot(A_minus(w_l), w_l)

def T_rot(beta):
    # this is the 3D rotation matrix (2.2.51 in Pavel's master thesis)
    # in 2D without the middle column and row, alpha = 0
    alpha = 0
    return array([
        [1, 0, 0, 0],
        [0, cos(alpha)*cos(beta), sin(beta), 0],
        [0, -cos(alpha)*sin(beta), cos(beta), 0],
        [0, 0, 0, 1]
        ])

w = array([1.1, -10, 13, 700.1])

print "R:\n", R(w)
print "R_inv:\n", R_inv(w)
print "D_minus:\n", D_minus(w)
#print "R_inv (numpy):\n", inv(R(w))
#print "R*R_inv:\n", dot(R(w), R_inv(w))
print "A:\n", A_minus(w)
print "-"*80
w_l = array([1.1, -10, 13, 700.1])
w_r = array([1.1, -10, 13, 800.1])
print "1st condition"
print "f_riemann(%r, %r):\n%r" % (w_l, w_l, f_riemann(w_l, w_l))
print f_x(w_l)
print "2nd condition"
alpha = 0.3
m = T_rot(alpha)
w_r = dot(m, w_r)
w_l = dot(m, w_l)
print "f_riemann 1:\n%r" % (f_riemann(w_r, w_l))
m = T_rot(pi)
print "f_riemann 2:\n%r" % (-dot(inv(m), f_riemann(dot(m, w_l), dot(m, w_r))))
