from __future__ import division

from pyparsing import (Word, Combine, Optional, alphas, alphanums, oneOf,
                delimitedList, Group, nums, Literal, OneOrMore,
                CaselessLiteral, Forward, ZeroOrMore, restOfLine,
                ParseException)

class ParseError(Exception):
    pass

# identifier
ident = Word( alphas, alphanums + '_' )

def expr_grammar():
    """
    Returns the grammar for an expression.
    """
    point = Literal( "." )
    e     = CaselessLiteral( "E" )
    fnumber = Combine( Word( "+-"+nums, nums ) +
                       Optional( point + Optional( Word( nums ) ) ) +
                       Optional( e + Word( "+-"+nums, nums ) ) )

    plus  = Literal( "+" )
    minus = Literal( "-" )
    mult  = Literal( "*" )
    div   = Literal( "/" )
    lpar  = Literal( "(" )
    rpar  = Literal( ")" )
    addop  = plus | minus
    multop = mult | div
    expop = Literal( "^" )
    pi    = CaselessLiteral( "pi" )

    expr = Forward()
    atom = (Optional("-") + ( pi | fnumber | ident + lpar + expr + rpar | ident ) | ( lpar + expr + rpar ))

    # by defining exponentiation as "atom [ ^ factor ]..." instead of "atom [ ^ atom ]...", we get right-to-left exponents, instead of left-to-righ
    # that is, 2^3^2 = 2^(3^2), not (2^3)^2.
    factor = Forward()
    factor << atom + ZeroOrMore( ( expop + factor ) )

    term = factor + ZeroOrMore( ( multop + factor ) )
    expr << term + ZeroOrMore( ( addop + term ) )
    return Combine(expr)

# now we define the grammar for the whole mesh file:
lbrace = Literal("{").suppress()
rbrace = Literal("}").suppress()
semicolon = Literal(";").suppress()
equal = Literal("=").suppress()
comment = Literal("#") + Optional(restOfLine)

expr = expr_grammar()
list_ = Forward()
item = expr | Group(list_)
list_ << (lbrace + delimitedList(item) + rbrace)
assig = ident + equal + item + Optional(semicolon)
mesh = OneOrMore(Group(assig))
mesh.ignore(comment)


def evaluate(s, namespace):
    """
    Evaluates the string "s" in the namespace of "namespace".

    The math module is automatically included in globals and "^" is converted
    to "**". Otherwise it has to be a valid Python syntax.
    """
    import math
    glob = math.__dict__
    s = s.replace("^", "**")
    try:
        r = eval(s, glob, namespace)
    except:
        raise ParseError("Failed to evaluate: %s" % s)
    return r

def evaluate_list(s, namespace):
    """
    Evaluates each item in the list recursively.

    Converts the list to a tuple.
    """
    if hasattr(s, "__iter__"):
        return tuple([evaluate_list(y, namespace) for y in s])
    else:
        if isinstance(s, str):
            return evaluate(s, namespace)
        else:
            return x

def read_hermes_format(filename):
    """
    Reads a mesh from a file in a hermes format.

    Returns nodes, elements, boundary, nurbs or raises a ParseError if the
    syntax is invalid.
    """
    m = open(filename).read()
    return read_hermes_format_str(m)

def read_hermes_format_str(m):
    """
    Reads a mesh from a string in a hermes format.

    Returns nodes, elements, boundary, nurbs or raises a ParseError if the
    syntax is invalid.
    """
    try:
        result = mesh.parseString(m)
    except ParseException, e:
        raise ParseError(str(e))
    namespace = {}
    for k, v in result:
        namespace[k] = evaluate_list(v, namespace)
    nodes = namespace.pop("vertices", None)
    elements = namespace.pop("elements", None)
    boundary = namespace.pop("boundaries", None)
    nurbs = namespace.pop("curves", None)
    if nodes is None or elements is None or boundary is None:
        raise ParseError("Either nodes, elements or boundary is missing")
    return nodes, elements, boundary, nurbs
