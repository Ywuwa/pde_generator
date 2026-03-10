# discretizer.py
from ast_nodes import Const, Var, Add, Mul, Func, Derivative, IndexedVar

def discretize(expr, constants):
    return _gen_cpp(expr, constants)


def _gen_cpp(expr, constants, shift=""):
    """
    Генерирует C++ код (используется для явных уравнений с частной производной по времени)
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

    raise NotImplementedError(f"Unsupported discretization: {type(expr)}")
    
# ==========================================
# 1. DISCRETIZE → AST
# ==========================================

def discretize_ast(expr, constants):

    if isinstance(expr, Const):
        return expr

    if isinstance(expr, Var):
        if expr.name in constants:
            return expr
        return IndexedVar(expr.name, "index")
    
    if isinstance(expr, IndexedVar):
      # IndexedVar уже дискретизирован, просто возвращаем как есть
      return expr

    if isinstance(expr, Add):
        return Add(
            discretize_ast(expr.left, constants),
            discretize_ast(expr.right, constants)
        )

    if isinstance(expr, Mul):
        return Mul(
            discretize_ast(expr.left, constants),
            discretize_ast(expr.right, constants)
        )

    if isinstance(expr, Func):
        return Func(expr.name, discretize_ast(expr.arg, constants))

    if isinstance(expr, Derivative):

        axis = expr.axis
        offset = f"offset_{axis}"
        h = f"h_{axis}"

        # ----------------------------------
        # Вторые производные
        # ----------------------------------

        if isinstance(expr.expr, Derivative):

            axis2 = expr.expr.axis
            inner = expr.expr.expr

            o1 = f"offset_{axis}"
            o2 = f"offset_{axis2}"

            pp = _shift(inner, f"index + {o1} + {o2}", constants)
            pm = _shift(inner, f"index + {o1} - {o2}", constants)
            mp = _shift(inner, f"index - {o1} + {o2}", constants)
            mm = _shift(inner, f"index - {o1} - {o2}", constants)

            numerator = Add(
                Add(pp, Mul(Const(-1), pm)),
                Add(Mul(Const(-1), mp), mm)
            )

            return Mul(
                numerator,
                Const(f"1/(4*{h}*h_{axis2})")
            )

        # ----------------------------------
        # Первая производная
        # ----------------------------------

        plus = _shift(expr.expr, f"index + {offset}", constants)
        minus = _shift(expr.expr, f"index - {offset}", constants)
        
        plus = expand_products(plus)
        minus = expand_products(minus)
        
        numerator = Add(
            plus,
            Mul(Const(-1), minus)
        )
        
        return Mul(
            numerator,
            Const(f"1/(2*{h})")
        )

    raise NotImplementedError(f"Unsupported discretization: {type(expr)}")

def expand_products(expr):

    if isinstance(expr, Mul):

        if isinstance(expr.left, Add):
            return Add(
                expand_products(Mul(expr.left.left, expr.right)),
                expand_products(Mul(expr.left.right, expr.right))
            )

        if isinstance(expr.right, Add):
            return Add(
                expand_products(Mul(expr.left, expr.right.left)),
                expand_products(Mul(expr.left, expr.right.right))
            )

        return Mul(
            expand_products(expr.left),
            expand_products(expr.right)
        )

    if isinstance(expr, Add):
        return Add(
            expand_products(expr.left),
            expand_products(expr.right)
        )

    return expr
  
# ==========================================
# SHIFT (добавляет индекс)
# ==========================================
def _shift(expr, index, constants):

    if isinstance(expr, Var):

        if expr.name in constants:
            return expr

        return IndexedVar(expr.name, index)

    if isinstance(expr, Add):
        return Add(
            _shift(expr.left, index, constants),
            _shift(expr.right, index, constants)
        )

    if isinstance(expr, Mul):
        return Mul(
            _shift(expr.left, index, constants),
            _shift(expr.right, index, constants)
        )

    if isinstance(expr, Func):
        return Func(expr.name, _shift(expr.arg, index, constants))

    return expr


# ==========================================
# 2. GENERATE C++
# ==========================================

def generate_cpp(expr):

    if isinstance(expr, Const):
        return str(expr.value)

    if isinstance(expr, Var):
        return expr.name

    if isinstance(expr, IndexedVar):
        return f"{expr.name}[{expr.index}]"

    if isinstance(expr, Add):
        return f"({generate_cpp(expr.left)} + {generate_cpp(expr.right)})"

    if isinstance(expr, Mul):
        return f"({generate_cpp(expr.left)} * {generate_cpp(expr.right)})"

    if isinstance(expr, Func):
        return f"{expr.name}({generate_cpp(expr.arg)})"

    raise NotImplementedError("Unsupported node")