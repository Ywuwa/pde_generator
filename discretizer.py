# discretizer.py
from ast_nodes import Const, Var, Add, Mul, Func, Derivative, TimeDerivative, IndexedVar, Scheme   
# ==========================================
# 1. DISCRETIZE → AST
# ==========================================

def discretize_ast(expr, constants):
    if expr is None:
      raise ValueError("expr is None before discretization")

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
        scheme = expr.scheme
    
        offset = f"offset_{axis}"
        h = f"h_{axis}"
    
        # ----------------------------------
        # Вторая производная (оставляем central)
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
    
        # CENTRAL
        if scheme == Scheme.CENTRAL:
            plus = _shift(expr.expr, f"index + {offset}", constants)
            minus = _shift(expr.expr, f"index - {offset}", constants)
    
            numerator = Add(plus, Mul(Const(-1), minus))
    
            return Mul(
                numerator,
                Const(f"1/(2*{h})")
            )
    
        # FORWARD
        if scheme == Scheme.FORWARD:
            plus = _shift(expr.expr, f"index + {offset}", constants)
            center = _shift(expr.expr, "index", constants)
    
            numerator = Add(plus, Mul(Const(-1), center))
    
            return Mul(
                numerator,
                Const(f"1/{h}")
            )
    
        # BACKWARD
        if scheme == Scheme.BACKWARD:
            center = _shift(expr.expr, "index", constants)
            minus = _shift(expr.expr, f"index - {offset}", constants)
    
            numerator = Add(center, Mul(Const(-1), minus))
    
            return Mul(
                numerator,
                Const(f"1/{h}")
            )
        raise ValueError(f"Unknown scheme: {scheme}")
      
    if isinstance(expr, TimeDerivative):
        inner = expr.expr
    
        if not isinstance(inner, Var):
            raise ValueError("Dt must be applied to a variable")
    
        new = IndexedVar(inner.name , "index")
        old = IndexedVar(inner.name + "0", "index")
    
        return Mul(
            Const("1/tau"),
            Add(new, Mul(Const(-1), old))
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