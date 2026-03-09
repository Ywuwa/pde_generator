from ast_nodes import Const, Var, Add, Mul, Func, Derivative, AXES
# calculus
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


# AST utilities
def contains_var(expr, var):
    """
    Ищет вхождение переменной var
    """
    if isinstance(expr, Var):
        return expr.name == var

    if isinstance(expr, Add):
        return contains_var(expr.left, var) or contains_var(expr.right, var)

    if isinstance(expr, Mul):
        return contains_var(expr.left, var) or contains_var(expr.right, var)

    if isinstance(expr, Func):
        return contains_var(expr.arg, var)

    if isinstance(expr, Derivative):
        return contains_var(expr.expr, var)

    return False
  
def flatten_add(expr):
    """
    Приводит дерево сложения в плоский список слагаемых
    """
    if isinstance(expr, Add):
        return flatten_add(expr.left) + flatten_add(expr.right)
    return [expr]

def build_sum(terms):
    """
    Строит сумму
    """
    if not terms:
        return None
    expr = terms[0]
    for t in terms[1:]:
        expr = Add(expr, t)
    return expr



# equation transforms
def split_by_variable(expr, variable):
    """
    Разносит слагаемые по разные стороны от знака равенства,
    в зависимости от вхождения variable в слагаемое
    """
    terms = flatten_add(expr)
    left = []
    right = []

    for t in terms:
        if contains_var(t, variable):
            left.append(t)
        else:
            right.append(t)

    lhs = build_sum(left)
    rhs = build_sum(right)
    return lhs, rhs



# algebra
def distributive_expand(expr: str):
    """
    Раскрывает A*(B+C) -> A*B + A*C на уровне AST
    """
    
    if isinstance(expr, Add):
        return Add(
            distributive_expand(expr.left),
            distributive_expand(expr.right)
        )

    if isinstance(expr, Mul):
        left = distributive_expand(expr.left)
        right = distributive_expand(expr.right)

        # A*(B+C)
        if isinstance(right, Add):
            return Add(
                distributive_expand(Mul(left, right.left)),
                distributive_expand(Mul(left, right.right))
            )

        # (A+B)*C
        if isinstance(left, Add):
            return Add(
                distributive_expand(Mul(left.left, right)),
                distributive_expand(Mul(left.right, right))
            )

        return Mul(left, right)

    return expr