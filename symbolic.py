from ast_nodes import Const, Var, Add, Mul, Func, Derivative, AXES, Scheme
# calculus
def diff(expr, axis, scheme=Scheme.CENTRAL, constants=None):
    """
    Дифференцирует AST-выражение по правилам
    """
    if isinstance(expr, Var):
      if constants and expr.name in constants:
          return Const(0)
      return Derivative(expr, axis, scheme)
    
    if isinstance(expr, Derivative):
        return Derivative(expr, axis, scheme)
    
    if isinstance(expr, Add):
        return Add(
            diff(expr.left, axis, scheme, constants),
            diff(expr.right, axis, scheme, constants)
        )
    
    if isinstance(expr, Mul):
        return Add(
            Mul(expr.left, diff(expr.right, axis, scheme, constants)),
            Mul(expr.right, diff(expr.left, axis, scheme, constants))
        )

    if isinstance(expr, Func):

        inner = diff(expr.arg, axis)

        if expr.name == "sin":
            return Mul(Func("cos", expr.arg), inner)

        if expr.name == "cos":
            return Mul(Mul(Const(-1), Func("sin", expr.arg)), inner)

        if expr.name == "exp":
            return Mul(Func("exp", expr.arg), inner)
          
        if expr.name == "log":
            return Mul(
                Mul(Const(1), 
                    Mul(Const(1), inner)),
                Mul(Const(1), 
                    Mul(Const(1), expr.arg))
            )

    raise NotImplementedError("Unsupported differentiation case")


def laplace(expr, scheme=Scheme.CENTRAL, constants=None):
    """
    Обработка лапласиана
    """
    result = None
    for axis in AXES:
        second = diff(diff(expr, axis, scheme, constants), axis, scheme, constants)
        result = second if result is None else Add(result, second)
    return result

def simplify_ast(expr):
    from ast_nodes import Const, Add, Mul

    if isinstance(expr, Add):
        l = simplify_ast(expr.left)
        r = simplify_ast(expr.right)

        # 0 + x
        if isinstance(l, Const) and l.value == 0:
            return r
        if isinstance(r, Const) and r.value == 0:
            return l

        return Add(l, r)

    if isinstance(expr, Mul):
        l = simplify_ast(expr.left)
        r = simplify_ast(expr.right)

        # 0 * x
        if (isinstance(l, Const) and l.value == 0) or \
           (isinstance(r, Const) and r.value == 0):
            return Const(0)

        # 1 * x
        if isinstance(l, Const) and l.value == 1:
            return r
        if isinstance(r, Const) and r.value == 1:
            return l

        return Mul(l, r)

    return expr