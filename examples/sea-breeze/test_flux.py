from numpy import array, zeros, dot
from numpy.linalg import inv

from _numerical_flux import matrix_R, matrix_R_inv

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

w = array([1.1, 10, 13, 700.1])

print "R:\n", R(w)
print "R_inv:\n", R_inv(w)
print "R_inv (numpy):\n", inv(R(w))
print "R*R_inv:\n", dot(R(w), R_inv(w))
