from numpy.linalg import solve

def run(sys):
    A = sys.get_matrix().todense()
    rhs = sys.get_rhs()
    print A
    print rhs
    print "-"*70
    print "solution:"
    print solve(A, rhs)
