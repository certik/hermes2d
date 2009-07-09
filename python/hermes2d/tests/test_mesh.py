from hermes2d.mesh import read_hermes_format_str, expr, ParseError
from hermes2d import raises

mesh1 = """\
vertices =
{
  { 0, -1 }, # first vertex
  { 1, -1 },
  { -1, 0 },
  { 0, 0 },
  { 1, 0 },
  { -1, 1 },
  { 0, 1 },
  { 0.707106781, 0.707106781 }
}

elements =
{
  { 0, 1, 4, 3, 0 },
  { 3, 4, 7, 0 },
  { 3, 7, 6, 0 },
  { 2, 3, 6, 5, 0 }
}

boundaries =
{
  { 0, 1, 1 },
  { 1, 4, 2 },
  { 3, 0, 4 },
  { 4, 7, 2 },
  { 7, 6, 2 },
  { 2, 3, 4 },
  { 6, 5, 2 },
  { 5, 2, 3 }
}

curves =
{
  { 4, 7, 45 },
  { 7, 6, 45 }
}
"""

mesh2 = """\
vertices =
{
  { 0.1, 0 },
  { 0.07071067809999999, 0.07071067809999999 },
  { 0, 0.1 },
  { 0.1707106781, 0 },
  { 0.1707106781, 0.07071067809999999 },
  { 0.1707106781, 0.1707106781 },
  { 0.07071067809999999, 0.1707106781 },
  { 0, 0.1707106781 },
  { 1, 0 },
  { 1, 0.07071067809999999 },
  { 1, 0.1707106781 },
  { 1, 1 },
  { 0.1707106781, 1 },
  { 0.07071067809999999, 1 },
  { 0, 1 },
  { -0.9, 1 },
  { -0.9, 0.1707106781 },
  { -0.9, 0.1 }
}

elements =
{
  { 4, 1, 0, 0 },
  { 6, 2, 1, 0 },
  { 3, 4, 0, 0 },
  { 6, 7, 2, 0 },
  { 1, 4, 5, 6, 0 },
  { 3, 8, 9, 4, 0 },
  { 4, 9, 10, 5, 0 },
  { 5, 10, 11, 12, 0 },
  { 6, 5, 12, 13, 0 },
  { 7, 6, 13, 14, 0 },
  { 16, 7, 14, 15, 0 },
  { 17, 2, 7, 16, 0 }
}

boundaries =
{
  { 1, 0, 5 },
  { 2, 1, 5 },
  { 0, 3, 1 },
  { 3, 8, 1 },
  { 8, 9, 2 },
  { 9, 10, 2 },
  { 10, 11, 2 },
  { 11, 12, 3 },
  { 12, 13, 3 },
  { 13, 14, 3 },
  { 14, 15, 3 },
  { 15, 16, 4 },
  { 17, 2, 5 },
  { 16, 17, 4 }
}

curves =
{
  { 1, 0, -45 },
  { 2, 1, -45 }
}
"""

mesh3 = """\
a = 1.0  # size of the mesh
b = sqrt(2)/2

vertices =
{
  { 0, -a },    # vertex 0
  { a, -a },    # vertex 1
  { -a, 0 },    # vertex 2
  { 0, 0 },     # vertex 3
  { a, 0 },     # vertex 4
  { -a, a },    # vertex 5
  { 0, a },     # vertex 6
  { a*b, a*b }  # vertex 7
}

elements =
{
  { 0, 1, 4, 3, 0 },  # quad 0
  { 3, 4, 7, 0 },     # tri 1
  { 3, 7, 6, 0 },     # tri 2
  { 2, 3, 6, 5, 0 }   # quad 3
}

boundaries =
{
  { 0, 1, 1 },
  { 1, 4, 2 },
  { 3, 0, 4 },
  { 4, 7, 2 },
  { 7, 6, 2 },
  { 2, 3, 4 },
  { 6, 5, 2 },
  { 5, 2, 3 }
}

curves =
{
  { 4, 7, 45 },  # +45 degree circular arcs
  { 7, 6, 45 }
}
"""

mesh4 = """\
a = 1.0  # horizontal size of the square
b = 1.0  # vertical size of the square

vertices =
{
  { 0, 0 },     # vertex 0
  { a, 0 },    # vertex 1
  { a, b },    # vertex 2
  { 0, b }     # vertex 3
}

elements =
{
  { 0, 1, 2, 3, 0 }  # quad 0
}

boundaries =
{
  { 0, 1, 1 },
  { 1, 2, 1 },
  { 2, 3, 1 },
  { 3, 0, 1 }
}
"""

mesh5 = """\
L = 15            # domain length (should be a multiple of 3)
H = 5             # domain height
S1 = 5/2          # x-center of circle
S2 = 5/2          # y-center of circle
R = 1             # circle radius
A = 1/(2*sqrt(2)) # helper length
EPS = 0.20        # inducing mesh non-symmetry

vertices =
{
  { 0, 0 },            # 0
  { S1 - A, 0 },       # 1
  { S1 + A, 0 },       # 2
  { L/3, 0 },          # 3
  { 2*L/3, 0 },        # 4
  { L, 0 },            # 5
  { 0, S2 - A },       # 6
  { S1 - A, S2 - A },  # 7
  { S1 + A, S2 - A },  # 8
  { L/3, S2 - A },     # 9
  { 2*L/3, S2 - A },   # 10
  { L, S2 - A },       # 11
  { 0, S2 + A },       # 12
  { S1 - A, S2 + A },  # 13
  { S1 + A, S2 + A },  # 14
  { L/3, S2 + A + EPS},  # 15 (inducing mesh non-symmetry)
  { 2*L/3, S2 + A },   # 16
  { L, S2 + A },       # 17
  { 0, H },            # 18
  { S1 - A, H },       # 19
  { S1 + A, H },       # 20
  { L/3, H },          # 21
  { 2*L/3, H },        # 22
  { L, H }             # 23
}

elements =
{
  { 0, 1, 7, 6, 0 },
  { 1, 2, 8, 7, 0 },
  { 2, 3, 9, 8, 0 },
  { 3, 4, 10, 9, 0 },
  { 4, 5, 11, 10, 0 },
  { 6, 7, 13, 12, 0 },
  { 8, 9, 15, 14, 0 },
  { 9, 10, 16, 15, 0 },
  { 10, 11, 17, 16, 0 },
  { 12, 13, 19, 18, 0 },
  { 13, 14, 20, 19, 0 },
  { 14, 15, 21, 20, 0 },
  { 15, 16, 22, 21, 0 },
  { 16, 17, 23, 22, 0 }
}

boundaries =
{
  { 0, 1, 1 },
  { 1, 2, 1 },
  { 2, 3, 1 },
  { 3, 4, 1 },
  { 4, 5, 1 },
  { 5, 11, 2 },
  { 11, 17, 2 },
  { 17, 23, 2 },
  { 23, 22, 3 },
  { 22, 21, 3 },
  { 21, 20, 3 },
  { 20, 19, 3 },
  { 19, 18, 3},
  { 18, 12, 4},
  { 12, 6, 4 },
  { 6, 0, 4 },
  { 7, 13, 5 },
  { 13, 14, 5},
  { 14, 8, 5 },
  { 8, 7, 5}
}

curves =
{
  { 7, 8, 90 },   # 45 degrees circular arc
  { 8, 14, 90 },  # 45 degrees circular arc
  { 14, 13, 90 }, # 45 degrees circular arc
  { 13, 7, 90 }   # 45 degrees circular arc
}

refinements =
{
  { 3, 2},
  { 0, 0},
  { 1, 1},
  { 2, 0},
  { 4, 2},
  { 5, 2},
  { 6, 2},
  { 7, 2},
  { 8, 2},
  { 9, 0},
  { 10, 1},
  { 11, 0},
  { 12, 2},
  { 13, 2},
  { 14, 0},
  { 15, 0},
  { 26, 0},
  { 27, 0},
  { 32, 2},
  { 33, 2},
  { 34, 2},
  { 35, 2},
  { 46, 0},
  { 47, 0},
  { 48, 0},
  { 49, 0}
}
"""

mesh6 = """\
a = 0.25   # horizontal size of an eleemnt
b = 0.1    # vertical size of an element
w = 0.000  # width of the cracks

vertices =
{
  { 0*a, 0 },         # vertex 0
  { 1*a, 0 },         # vertex 1
  { 2*a, 0 },         # vertex 2
  { 3*a, 0 },         # vertex 3
  { 4*a, 0 },         # vertex 4
  { 5*a, 0 },         # vertex 5
  { 6*a, 0 },         # vertex 6
  { 0*a, b },         # vertex 7
  { 1*a, b },         # vertex 8
  { 2*a, b },         # vertex 9
  { 3*a, b },         # vertex 10
  { 4*a, b - w/2},    # vertex 11
  { 4*a, b + w/2},    # vertex 12
  { 5*a, b },         # vertex 13
  { 6*a, b },         # vertex 14
  { 0*a, 2*b },       # vertex 15
  { 1*a, 2*b },       # vertex 16
  { 2*a, 2*b - w/2},  # vertex 17
  { 2*a, 2*b + w/2},  # vertex 18
  { 3*a, 2*b },       # vertex 19
  { 4*a, 2*b },       # vertex 20
  { 5*a, 2*b },       # vertex 21
  { 6*a, 2*b },       # vertex 22
  { 0*a, 3*b },       # vertex 23
  { 1*a, 3*b },       # vertex 24
  { 2*a, 3*b },       # vertex 25
  { 3*a, 3*b },       # vertex 26
  { 4*a, 3*b },       # vertex 27
  { 5*a, 3*b },       # vertex 28
  { 6*a, 3*b }        # vertex 29
}

elements =
{
  { 0, 1, 8, 7, 0 },      # quad 0
  { 1, 2, 9, 8, 0 },      # quad 1
  { 2, 3, 10, 9, 0 },     # quad 2
  { 3, 4, 11, 10, 0 },    # quad 3
  { 4, 5, 13, 11, 0 },    # quad 4
  { 5, 6, 14, 13, 0 },    # quad 5

  { 7, 8, 16, 15, 0 },    # quad 6
  { 8, 9, 17, 16, 0 },    # quad 7
  { 9, 10, 19, 17, 0 },   # quad 8
  { 10, 12, 20, 19, 0 },  # quad 9
  { 12, 13, 21, 20, 0 },  # quad 10
  { 13, 14, 22, 21, 0 },  # quad 11

  { 15, 16, 24, 23, 0 },  # quad 12
  { 16, 18, 25, 24, 0 },  # quad 13
  { 18, 19, 26, 25, 0 },  # quad 14
  { 19, 20, 27, 26, 0 },  # quad 15
  { 20, 21, 28, 27, 0 },  # quad 16
  { 21, 22, 29, 28, 0 }   # quad 17
}

boundaries =
{
  { 0, 1, 3 },   # bottom edges
  { 1, 2, 3 },
  { 2, 3, 3 },
  { 3, 4, 3 },
  { 4, 5, 3 },
  { 5, 6, 3 },

  { 6, 14, 3 },  # right edges
  { 14, 22, 3 },
  { 22, 29, 3 },

  { 29, 28, 2 }, # top edges
  { 28, 27, 2 },
  { 27, 26, 2 },
  { 26, 25, 2 },
  { 25, 24, 2 },
  { 24, 23, 2 },


  { 23, 15, 1 }, # left edges
  { 15, 7, 1 },
  { 7, 0, 1 },

  { 16, 17, 3 }, # left crack
  { 17, 19, 3 },
  { 19, 18, 3 },
  { 18, 16, 3 },

  { 10, 11, 3 }, # right crack
  { 11, 13, 3 },
  { 13, 12, 3 },
  { 12, 10, 3 }


}
"""

mesh7 = """ \
t = 0.1  # thickness
l = 0.7  # length

left = 1;
top  = 2;
rest = 3;


a = sqrt(l^2 - (l-t)^2)
b = t
alpha = atan(b/l)
delta = atan(a/(l-t))
beta  = delta - alpha
gamma = pi/2 - 2*delta
c = (l-t)*sin(alpha)
d = (l-t)*cos(alpha)
e = (l-t)*sin(delta)
f = (l-t)*cos(delta)
q = sqrt(2)/2


vertices =
{
  { l-t, 0 },  # 0
  { l, 0 },    # 1
  { d, c },    # 2
  { l, b },    # 3
  { f, e },    # 4
  { l-t, a },  # 5
  { l, a },    # 6

  { 0, l-t },  # 7
  { 0, l },    # 8
  { c, d },    # 9
  { b, l },    # 10
  { e, f },    # 11
  { a, l-t },  # 12
  { a, l },    # 13

  { l-t, l-t }, # 14
  { l, l-t },   # 15
  { l, l },     # 16
  { l-t, l },   # 17

  { l, -t },       # 18
  { l-q*t, -q*t }, # 19
  { -t, l },       # 20
  { -q*t, l-q*t }  # 21
}


m = 0

elements =
{
  { 0, 1, 3, 2, m },
  { 2, 3, 5, 4, m },
  { 6, 5, 3, m },
  { 8, 7, 9, 10, m },
  { 10, 9, 11, 12, m },
  { 13, 10, 12, m },
  { 4, 5, 12, 11, m },
  { 5, 6, 15, 14, m },
  { 13, 12, 14, 17, m },
  { 14, 15, 16, 17, m },
  { 0, 19, 1, m },
  { 19, 18, 1, m },
  { 21, 7, 8, m },
  { 20, 21, 8, m }
}

boundaries =
{
  { 18, 1, left },
  { 1, 3, left },
  { 3, 6, left },
  { 6, 15, left },
  { 15, 16, left },
  { 16, 17, top },
  { 17, 13, top },
  { 13, 10, top },
  { 10, 8, top },
  { 8, 20, top },
  { 20, 21, rest },
  { 21, 7, rest },
  { 7, 9, rest },
  { 9, 11, rest },
  { 11, 4, rest },
  { 4, 2, rest },
  { 2, 0, rest },
  { 0, 19, rest },
  { 19, 18, rest },
  { 5, 14, rest },
  { 14, 12, rest },
  { 12, 5, rest }
}


alpha = 180*alpha/pi
beta  = 180*beta/pi
gamma = 180*gamma/pi

curves =
{
  { 0, 2, alpha },
  { 2, 4, beta },
  { 4, 11, gamma },
  { 11, 9, beta },
  { 9, 7, alpha },
  { 5,12, gamma },
  { 0, 19, 45.0 },
  { 19, 18, 45.0 },
  { 20, 21, 45.0 },
  { 21, 7, 45.0 }
};
"""

def flat(x):
    if hasattr(x, "__iter__"):
        r = []
        for y in x:
            r.extend(flat(y))
        return r
    else:
        return [x]

def compare(a, b, eps=1e-9):
    a = flat(a)
    b = flat(b)
    for x, y in zip(a, b):
        if not (abs(x-y) < eps):
            return False
    return True

def test_expr():
    def t(x, correct_result=None):
        r = list(expr.parseString(x))
        if correct_result is not None:
            return r == [correct_result]
        else:
            return r == [x]
    assert t("1")
    assert t("-1")
    assert t("+1")
    assert t("1.245")
    assert t("-1.245")
    assert t("+1.245")
    assert t("a")
    assert t("-a")
    assert t("a-b")
    assert t("sqrt(2)/2")
    assert t("a*b")
    assert t("5/2")
    assert t("1/(2*sqrt(2))")
    assert t("0.20")
    assert t("S1 - A", "S1-A")
    assert t("S1 + A", "S1+A")
    assert t("L/3")
    assert t("2*L/3")
    assert t("L")
    assert t("S2 - A", "S2-A")
    assert t("S2 + A + EPS", "S2+A+EPS")
    assert t("4*a")
    assert t("b - w/2", "b-w/2")
    assert t("b + w/2", "b+w/2")
    assert t("0.1")
    assert t("0.7")
    assert t("sqrt(l^2 - (l-t)^2)", "sqrt(l^2-(l-t)^2)")
    assert t("atan(b/l)")
    assert t("atan(a/(l-t))")
    assert t("delta - alpha", "delta-alpha")
    assert t("pi/2 - 2*delta", "pi/2-2*delta")
    assert t("(l-t)*sin(alpha)")
    assert t("(l-t)*cos(alpha)")
    assert t("(l-t)*sin(delta)")
    assert t("(l-t)*cos(delta)")
    assert t("180*alpha/pi")
    assert t("180*beta/pi")
    assert t("180*gamma/pi")
    assert t("l-t")
    assert t("l-q*t")
    assert t("-q*t")
    assert t("A*sin(x) + B*cos(x) + pi/2 - some_result",
            "A*sin(x)+B*cos(x)+pi/2-some_result")

def test_loader1():
    nodes, elements, boundary, nurbs = read_hermes_format_str(mesh1)
    assert compare(nodes, [[0.0, -1.0], [1.0, -1.0], [-1.0, 0.0], [0.0, 0.0],
        [1.0, 0.0], [-1.0, 1.0], [0.0, 1.0], [0.70710678100000002,
            0.70710678100000002]])
    assert elements == ( ( 0, 1, 4, 3, 0 ), ( 3, 4, 7, 0 ),
        ( 3, 7, 6, 0 ),
        ( 2, 3, 6, 5, 0 ))
    assert boundary == (
          ( 0, 1, 1 ),
          ( 1, 4, 2 ),
          ( 3, 0, 4 ),
          ( 4, 7, 2 ),
          ( 7, 6, 2 ),
          ( 2, 3, 4 ),
          ( 6, 5, 2 ),
          ( 5, 2, 3 )
          )
    assert nurbs == (
                  ( 4, 7, 45 ),
                  ( 7, 6, 45 )
                )
    assert boundary[1][0] == 1
    assert boundary[1][1] == 4
    assert boundary[1][2] == 2
    # this tests that ints are converted to ints, but not float
    assert isinstance(boundary[1][0], (int, long))
    assert isinstance(boundary[1][1], (int, long))
    assert isinstance(boundary[1][2], (int, long))

def test_loader2():
    nodes, elements, boundary, nurbs = read_hermes_format_str(mesh2)
    assert compare(nodes, [
          [ 0.1, 0 ],
          [ 0.07071067809999999, 0.07071067809999999 ],
          [ 0, 0.1 ],
          [ 0.1707106781, 0 ],
          [ 0.1707106781, 0.07071067809999999 ],
          [ 0.1707106781, 0.1707106781 ],
          [ 0.07071067809999999, 0.1707106781 ],
          [ 0, 0.1707106781 ],
          [ 1, 0 ],
          [ 1, 0.07071067809999999 ],
          [ 1, 0.1707106781 ],
          [ 1, 1 ],
          [ 0.1707106781, 1 ],
          [ 0.07071067809999999, 1 ],
          [ 0, 1 ],
          [ -0.9, 1 ],
          [ -0.9, 0.1707106781 ],
          [ -0.9, 0.1 ]
            ] )
    assert elements == (
          ( 4, 1, 0, 0 ),
          ( 6, 2, 1, 0 ),
          ( 3, 4, 0, 0 ),
          ( 6, 7, 2, 0 ),
          ( 1, 4, 5, 6, 0 ),
          ( 3, 8, 9, 4, 0 ),
          ( 4, 9, 10, 5, 0 ),
          ( 5, 10, 11, 12, 0 ),
          ( 6, 5, 12, 13, 0 ),
          ( 7, 6, 13, 14, 0 ),
          ( 16, 7, 14, 15, 0 ),
          ( 17, 2, 7, 16, 0 )
            )
    assert boundary == (
          ( 1, 0, 5 ),
          ( 2, 1, 5 ),
          ( 0, 3, 1 ),
          ( 3, 8, 1 ),
          ( 8, 9, 2 ),
          ( 9, 10, 2 ),
          ( 10, 11, 2 ),
          ( 11, 12, 3 ),
          ( 12, 13, 3 ),
          ( 13, 14, 3 ),
          ( 14, 15, 3 ),
          ( 15, 16, 4 ),
          ( 17, 2, 5 ),
          ( 16, 17, 4 )
          )
    assert nurbs == (
                  ( 1, 0, -45 ),
                  ( 2, 1, -45 )
                )

def test_loader3():
    nodes, elements, boundary, nurbs = read_hermes_format_str(mesh3)
    assert compare(nodes, [[0.0, -1.0], [1.0, -1.0], [-1.0, 0.0], [0.0, 0.0],
        [1.0, 0.0], [-1.0, 1.0], [0.0, 1.0], [0.70710678100000002,
            0.70710678100000002]])
    assert elements == ( ( 0, 1, 4, 3, 0 ), ( 3, 4, 7, 0 ),
        ( 3, 7, 6, 0 ),
        ( 2, 3, 6, 5, 0 ))
    assert boundary == (
          ( 0, 1, 1 ),
          ( 1, 4, 2 ),
          ( 3, 0, 4 ),
          ( 4, 7, 2 ),
          ( 7, 6, 2 ),
          ( 2, 3, 4 ),
          ( 6, 5, 2 ),
          ( 5, 2, 3 )
          )
    assert nurbs == (
                  ( 4, 7, 45 ),
                  ( 7, 6, 45 )
                )
    assert boundary[1][0] == 1
    assert boundary[1][1] == 4
    assert boundary[1][2] == 2
    # this tests that ints are converted to ints, but not float
    assert isinstance(boundary[1][0], (int, long))
    assert isinstance(boundary[1][1], (int, long))
    assert isinstance(boundary[1][2], (int, long))

def test_loader4():
    nodes, elements, boundary, nurbs = read_hermes_format_str(mesh4)
    assert compare(nodes, ((0, 0), (1.0, 0), (1.0, 1.0), (0, 1.0)))
    assert elements == ((0, 1, 2, 3, 0),)
    assert boundary == ((0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1))
    assert nurbs is None

def test_loader5():
    nodes, elements, boundary, nurbs = read_hermes_format_str(mesh5)
    assert compare(nodes, ((0, 0), (2.146446609406726, 0),
        (2.853553390593274, 0), (5.0, 0), (10.0, 0), (15, 0),
        (0, 2.146446609406726), (2.146446609406726, 2.146446609406726),
        (2.853553390593274, 2.146446609406726), (5.0, 2.146446609406726),
        (10.0, 2.146446609406726), (15, 2.146446609406726),
        (0, 2.853553390593274), (2.146446609406726, 2.853553390593274),
        (2.853553390593274, 2.853553390593274), (5.0, 3.0535533905932741),
        (10.0, 2.853553390593274), (15, 2.853553390593274), (0, 5),
        (2.146446609406726, 5), (2.853553390593274, 5), (5.0, 5), (10.0, 5),
        (15, 5)))
    assert elements == ((0, 1, 7, 6, 0), (1, 2, 8, 7, 0), (2, 3, 9, 8, 0),
            (3, 4, 10, 9, 0), (4, 5, 11, 10, 0), (6, 7, 13, 12, 0),
            (8, 9, 15, 14, 0), (9, 10, 16, 15, 0), (10, 11, 17, 16, 0),
            (12, 13, 19, 18, 0), (13, 14, 20, 19, 0), (14, 15, 21, 20, 0),
            (15, 16, 22, 21, 0), (16, 17, 23, 22, 0))
    assert boundary == ((0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1),
            (5, 11, 2), (11, 17, 2), (17, 23, 2), (23, 22, 3), (22, 21, 3),
            (21, 20, 3), (20, 19, 3), (19, 18, 3), (18, 12, 4), (12, 6, 4),
            (6, 0, 4), (7, 13, 5), (13, 14, 5), (14, 8, 5), (8, 7, 5))
    assert nurbs == ((7, 8, 90), (8, 14, 90), (14, 13, 90), (13, 7, 90))

def test_loader6():
    nodes, elements, boundary, nurbs = read_hermes_format_str(mesh6)
    assert compare(nodes, ((0.0, 0), (0.25, 0), (0.5, 0), (0.75, 0), (1.0, 0),
        (1.25, 0), (1.5, 0), (0.0, 0.1), (0.25, 0.1), (0.5, 0.1), (0.75, 0.1),
        (1.0, 0.1), (1.0, 0.1), (1.25, 0.1), (1.5, 0.1), (0.0, 0.2),
        (0.25, 0.2), (0.5, 0.2), (0.5, 0.2), (0.75, 0.2), (1.0, 0.2),
        (1.25, 0.2), (1.5, 0.2), (0.0, 0.3), (0.25, 0.3), (0.5, 0.3),
        (0.75, 0.3), (1.0, 0.3), (1.25, 0.3), (1.5, 0.3)))
    assert elements == ((0, 1, 8, 7, 0), (1, 2, 9, 8, 0), (2, 3, 10, 9, 0),
            (3, 4, 11, 10, 0), (4, 5, 13, 11, 0), (5, 6, 14, 13, 0),
            (7, 8, 16, 15, 0), (8, 9, 17, 16, 0), (9, 10, 19, 17, 0),
            (10, 12, 20, 19, 0), (12, 13, 21, 20, 0), (13, 14, 22, 21, 0),
            (15, 16, 24, 23, 0), (16, 18, 25, 24, 0), (18, 19, 26, 25, 0),
            (19, 20, 27, 26, 0), (20, 21, 28, 27, 0), (21, 22, 29, 28, 0))
    assert boundary == ((0, 1, 3), (1, 2, 3), (2, 3, 3), (3, 4, 3), (4, 5, 3),
            (5, 6, 3), (6, 14, 3), (14, 22, 3), (22, 29, 3), (29, 28, 2),
            (28, 27, 2), (27, 26, 2), (26, 25, 2), (25, 24, 2), (24, 23, 2),
            (23, 15, 1), (15, 7, 1), (7, 0, 1), (16, 17, 3), (17, 19, 3),
            (19, 18, 3), (18, 16, 3), (10, 11, 3), (11, 13, 3), (13, 12, 3),
            (12, 10, 3))
    assert nurbs is None

def test_loader7():
    nodes, elements, boundary, nurbs = read_hermes_format_str(mesh7)
    assert compare(nodes, ((0.6, 0), (0.7, 0),
        (0.59396969619669993, 0.084852813742385721), (0.7, 0.1),
        (0.51428571428571435, 0.30904725218262763), (0.6, 0.36055512754639885),
        (0.7, 0.36055512754639885), (0, 0.6), (0, 0.7),
        (0.084852813742385721, 0.59396969619669993), (0.1, 0.7),
        (0.30904725218262763, 0.51428571428571435), (0.36055512754639885, 0.6),
        (0.36055512754639885, 0.7), (0.6, 0.6), (0.7, 0.6), (0.7, 0.7),
        (0.6, 0.7), (0.7, -0.1), (0.62928932188134523, -0.070710678118654766),
        (-0.1, 0.7), (-0.070710678118654766, 0.62928932188134523)))
    assert elements == ((0, 1, 3, 2, 0), (2, 3, 5, 4, 0), (6, 5, 3, 0),
            (8, 7, 9, 10, 0), (10, 9, 11, 12, 0), (13, 10, 12, 0),
            (4, 5, 12, 11, 0), (5, 6, 15, 14, 0), (13, 12, 14, 17, 0),
            (14, 15, 16, 17, 0), (0, 19, 1, 0), (19, 18, 1, 0), (21, 7, 8, 0),
            (20, 21, 8, 0))
    assert boundary == ((18, 1, 1), (1, 3, 1), (3, 6, 1), (6, 15, 1),
            (15, 16, 1), (16, 17, 2), (17, 13, 2), (13, 10, 2), (10, 8, 2),
            (8, 20, 2), (20, 21, 3), (21, 7, 3), (7, 9, 3), (9, 11, 3),
            (11, 4, 3), (4, 2, 3), (2, 0, 3), (0, 19, 3), (19, 18, 3),
            (5, 14, 3), (14, 12, 3), (12, 5, 3))
    assert nurbs == ((0, 2, 8.13010235415598), (2, 4, 22.872616779718008),
            (4, 11, 27.994561732252027), (11, 9, 22.872616779718008),
            (9, 7, 8.13010235415598), (5, 12, 27.994561732252027),
            (0, 19, 45.0), (19, 18, 45.0), (20, 21, 45.0), (21, 7, 45.0))

def test_errors():
    mesh = """\
    a = %s
    vertices =
    {
      { 0, -1 }, # first vertex
      { 0.707106781, 0.707106781 }
    }

    elements =
    {
      { 0, 1, 4, 3, 0 }
    }

    boundaries =
    {
      { 0, 1, 1 }
    }
    """
    # this works:
    read_hermes_format_str(mesh % "34")
    # this fails:
    assert raises(ParseError, "read_hermes_format_str(mesh % 'x x x')")
    assert raises(ParseError, "read_hermes_format_str(mesh % '3 3 3')")
    assert raises(ParseError, "read_hermes_format_str(mesh % '+')")
    assert raises(ParseError, "read_hermes_format_str(mesh % '^^')")

def test_division():
    """
    Tests true division.

    E.g. 3/2 is 1.5, not 1 in the mesh file. This test tests it.
    """
    mesh = """\
    vertices =
    {
      { 0, -1 },
      { %s, 0.707106781 }
    }

    elements =
    {
      { 0, 1, 4, 3, 0 }
    }

    boundaries =
    {
      { 0, 1, 1 }
    }
    """
    nodes, elements, boundaries, _ = read_hermes_format_str(mesh % "3")
    assert compare(nodes, ((0, -1), (3, 0.70710678100000002)))

    nodes, elements, boundaries, _ = read_hermes_format_str(mesh % "3/2")
    assert compare(nodes, ((0, -1), (1.5, 0.70710678100000002)))
