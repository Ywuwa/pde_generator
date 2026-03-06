from ast_nodes import Const, Var, Add, Mul, Func, Derivative, AXES

def diff(expr, axis):
    """
    Дифференцирует AST-выражение по правилам
    """
    if isinstance(expr, Const):
        return Const(0)

    if isinstance(expr, Var):
        return Derivative(expr, axis)

    if isinstance(expr, Derivative):
        return Derivative(expr, axis)

    if isinstance(expr, Add):
        return Add(diff(expr.left, axis),
                   diff(expr.right, axis))

    if isinstance(expr, Mul):
        return Add(
            Mul(expr.left, diff(expr.right, axis)),
            Mul(expr.right, diff(expr.left, axis))
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


def laplace(expr):
    """
    Обработка лапласиана
    """
    result = None
    for axis in AXES:
        second = diff(diff(expr, axis), axis)
        result = second if result is None else Add(result, second)
    return result