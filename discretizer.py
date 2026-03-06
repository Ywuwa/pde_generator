# discretizer.py
from ast_nodes import Const, Var, Add, Mul, Func, Derivative

def discretize(expr, constants):
    return _gen_cpp(expr, constants)


def _gen_cpp(expr, constants, shift=""):
    """
    Генерирует C++ код
    """
    if isinstance(expr, Const):
        return str(expr.value)

    if isinstance(expr, Var):

        if expr.name in constants:
            return expr.name

        return f"{expr.name}[index{shift}]"

    if isinstance(expr, Add):
        return f"({_gen_cpp(expr.left, constants, shift)} + {_gen_cpp(expr.right, constants, shift)})"

    if isinstance(expr, Mul):
        return f"({_gen_cpp(expr.left, constants, shift)} * {_gen_cpp(expr.right, constants, shift)})"

    if isinstance(expr, Func):
        return f"{expr.name}({_gen_cpp(expr.arg, constants, shift)})"

    if isinstance(expr, Derivative):

        axis = expr.axis
        offset = f"offset_{axis}"
        h = f"h_{axis}"

        if isinstance(expr.expr, Derivative):

            axis2 = expr.expr.axis
            inner = expr.expr.expr

            o1 = f"offset_{axis}"
            o2 = f"offset_{axis2}"
            h1 = f"h_{axis}"
            h2 = f"h_{axis2}"

            pp = _gen_cpp(inner, constants, f" + {o1} + {o2}")
            pm = _gen_cpp(inner, constants, f" + {o1} - {o2}")
            mp = _gen_cpp(inner, constants, f" - {o1} + {o2}")
            mm = _gen_cpp(inner, constants, f" - {o1} - {o2}")

            return f"({pp} - {pm} - {mp} + {mm}) / (4*{h1}*{h2})"

        plus = _gen_cpp(expr.expr, constants, f" + {offset}")
        minus = _gen_cpp(expr.expr, constants, f" - {offset}")

        return f"({plus} - {minus}) / (2*{h})"

    raise NotImplementedError("Unsupported discretization")