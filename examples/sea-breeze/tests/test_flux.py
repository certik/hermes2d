from math import sin, cos, pi

from numpy import array, zeros, dot, eye
from numpy.linalg import inv

from _numerical_flux import matrix_R, matrix_R_inv, matrix_D_minus, \
        c_v, flux_riemann
from _numerical_flux import R as R_const

eps = 1e-10

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

def test_inv():
    w = array([1.1, -10, 13, 700.1])
    assert (dot(R(w), R_inv(w))-eye(4) < eps).all()
    assert (R_inv(w) - inv(R(w)) < eps).all()
    w = array([1.1, 10, 13, 700.1])
    assert (dot(R(w), R_inv(w))-eye(4) < eps).all()
    assert (R_inv(w) - inv(R(w)) < eps).all()
    w = array([3.1, 10, 13, 800.1])
    assert (dot(R(w), R_inv(w))-eye(4) < eps).all()
    assert (R_inv(w) - inv(R(w)) < eps).all()

def test_flux():
    w_l = array([1.1, -10, 13, 700.1])
    w_r = array([1.1, -10, 13, 800.1])
    assert (f_riemann(w_l, w_l) - f_x(w_l) < eps).all()
    assert (f_riemann(w_r, w_r) - f_x(w_r) < eps).all()

    alpha = 0
    m = T_rot(alpha)
    w_r = dot(m, w_r)
    w_l = dot(m, w_l)
    m = T_rot(pi)
    flux1 = f_riemann(w_l, w_r)
    flux2 = -dot(inv(m), f_riemann(dot(m, w_r), dot(m, w_l)))
    flux3 = flux_riemann(w_l, w_r)
    flux4 = -dot(inv(m), flux_riemann(dot(m, w_r), dot(m, w_l)))
    assert (flux1 - flux2 < eps).all()
    assert (flux1 - flux3 < eps).all()
    assert (flux1 - flux4 < eps).all()

    w_l = array([1.1, -10, 13, 700.1])
    w_r = array([1.1, -10, 13, 800.1])
    alpha = 0.3
    m = T_rot(alpha)
    w_r = dot(m, w_r)
    w_l = dot(m, w_l)
    m = T_rot(pi)
    flux1 = f_riemann(w_l, w_r)
    flux2 = -dot(inv(m), f_riemann(dot(m, w_r), dot(m, w_l)))
    flux3 = flux_riemann(w_l, w_r)
    flux4 = -dot(inv(m), flux_riemann(dot(m, w_r), dot(m, w_l)))
    assert (flux1 - flux2 < eps).all()
    assert (flux1 - flux3 < eps).all()
    assert (flux1 - flux4 < eps).all()
