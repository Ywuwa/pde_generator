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
            inner_der = expr.expr
        
            axis1 = expr.axis
            axis2 = inner_der.axis
        
            scheme1 = expr.scheme
            scheme2 = inner_der.scheme
        
            inner = inner_der.expr
        
            h1 = f"h_{axis1}"
            h2 = f"h_{axis2}"
        
            o1 = f"offset_{axis1}"
            o2 = f"offset_{axis2}"
            # Central-Central
            if scheme1 == Scheme.CENTRAL and scheme2 == scheme1:
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
                  Const(f"1/(4*{h1}*{h2})")
              )
            
            # Forward-Forward
            if scheme1 == Scheme.FORWARD and scheme2 == scheme1:
              pp = _shift(inner, f"index + {o1} + {o2}", constants)
              p1 = _shift(inner, f"index + {o1}", constants)
              p2 = _shift(inner, f"index + {o2}", constants)
              c = _shift(inner, "index", constants)
          
              numerator = Add(
                  Add(pp, c),
                  Mul(Const(-1), Add(p1, p2))
              )
          
              return Mul(
                  numerator,
                  Const(f"1/({h1}*{h2})")
              )
            
            # Backward-Backward
            if scheme1 == Scheme.BACKWARD and scheme2 == scheme1:
              mm = _shift(inner, f"index - {o1} - {o2}", constants)
              m1 = _shift(inner, f"index - {o1}", constants)
              m2 = _shift(inner, f"index - {o2}", constants)
              c = _shift(inner, "index", constants)
          
              numerator = Add(
                  Add(c, mm),
                  Mul(Const(-1), Add(m1, m2))
              )
          
              return Mul(
                  numerator,
                  Const(f"1/({h1}*{h2})")
              )
            
            # Forward-Backward
            if scheme1 == Scheme.FORWARD and scheme2 == Scheme.BACKWARD:
              pm = _shift(inner, f"index + {o1} - {o2}", constants)
              p = _shift(inner, f"index + {o1}", constants)
              m = _shift(inner, f"index - {o2}", constants)
              c = _shift(inner, "index", constants)
          
              numerator = Add(
                  Add(p, m),
                  Mul(Const(-1), Add(pm, c))
              )
          
              return Mul(
                  numerator,
                  Const(f"1/({h1}*{h2})")
              )
            
            # Backward-Forward
            if scheme1 == Scheme.BACKWARD and scheme2 == Scheme.FORWARD:
              mp = _shift(inner, f"index - {o1} + {o2}", constants)
              m = _shift(inner, f"index - {o1}", constants)
              p = _shift(inner, f"index + {o2}", constants)
              c = _shift(inner, "index", constants)
          
              numerator = Add(
                  Add(p, m),
                  Mul(Const(-1), Add(c, mp))
              )    
              return Mul(
                  numerator,
                  Const(f"1/({h1}*{h2})")
              )
            
            # Central-Forward
            if scheme1 == Scheme.CENTRAL and scheme2 == Scheme.FORWARD:
              pp = _shift(inner, f"index + {o1} + {o2}", constants)
              mp = _shift(inner, f"index - {o1} + {o2}", constants)
              p = _shift(inner, f"index + {o1}", constants)
              m = _shift(inner, f"index - {o1}", constants)
              
              numerator = Add(
                  Add(pp, m),
                  Mul(Const(-1), Add(mp, p))
              )
              return Mul(
                  numerator,
                  Const(f"1/(2*{h1}*{h2})")
              )
            
            # Forward-Central
            if scheme1 == Scheme.FORWARD and scheme2 == Scheme.CENTRAL:
              pp = _shift(inner, f"index + {o1} + {o2}", constants)
              pm = _shift(inner, f"index + {o1} - {o2}", constants)
              cp = _shift(inner, f"index + {o2}", constants)
              cm = _shift(inner, f"index - {o2}", constants)
              
              numerator = Add(
                  Add(pp, cm),
                  Mul(Const(-1), Add(pm, cp))
              )
              return Mul(
                  numerator,
                  Const(f"1/(2*{h1}*{h2})")
              )
            
            # Central-Backward
            if scheme1 == Scheme.CENTRAL and scheme2 == Scheme.BACKWARD:
              mm = _shift(inner, f"index - {o1} - {o2}", constants)
              pc = _shift(inner, f"index + {o1}", constants)
              mc = _shift(inner, f"index - {o1}", constants)
              pm = _shift(inner, f"index + {o1} - {o2}", constants)
              
              numerator = Add(
                  Add(mm, pc),
                  Mul(Const(-1), Add(pm, mc))
              )
              return Mul(
                  numerator,
                  Const(f"1/(2*{h1}*{h2})")
              )
            # Backward-Central
            if scheme1 == Scheme.BACKWARD and scheme2 == Scheme.CENTRAL:
              mm = _shift(inner, f"index - {o1} - {o2}", constants)
              cp = _shift(inner, f"index + {o1}", constants)
              cm = _shift(inner, f"index - {o1}", constants)
              mp = _shift(inner, f"index + {o1} - {o2}", constants)
              
              numerator = Add(
                  Add(mm, cp),
                  Mul(Const(-1), Add(mp, cm))
              )
              return Mul(
                  numerator,
                  Const(f"1/(2*{h1}*{h2})")
              )
            raise NotImplementedError(
                f"Unsupported second derivative scheme: {scheme1}, {scheme2}"
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
    
def simplify(expr):
  """
  expr : Expr
    AST-выражение. 
    Возвращает упрощённое AST-выражение: убирает нулевые константы, слагаемые и т.д.
  """
  if isinstance(expr, Add):
      left = simplify(expr.left)
      right = simplify(expr.right)
  
      # одинаковые константы
      if isinstance(left, Const) and isinstance(right, Const):
          return Const(left.value + right.value)
  
      # x + (-x)
      if str(left) == str(Mul(Const(-1), right)):
          return Const(0)
  
      if str(right) == str(Mul(Const(-1), left)):
          return Const(0)
  
      return Add(left, right)

  if isinstance(expr, Mul):
      left = simplify(expr.left)
      right = simplify(expr.right)
  
      if isinstance(left, Const) and isinstance(right, Const):
          return Const(left.value * right.value)
  
      return Mul(left, right)

  return expr