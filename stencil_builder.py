# stencil_builder.py
from ast_nodes import Const, Var, Add, Mul, IndexedVar, Derivative, TimeDerivative, Scheme
from stencil import StencilExpr
from symbolic import diff

def to_stencil(expr, constants):
    """
    Возвращает:
        { var_name: StencilExpr }
    """
    #print(expr)
    # --- Числовая константа ---
    if isinstance(expr, Const):
        return {"__const__": make_const_stencil(expr)}

    # --- Переменная ---
    if isinstance(expr, Var):
        # если это константа → используем имя, а не значение
        if expr.name in constants:
            return {
                "__const__": make_const_stencil(Const(expr.name))
            }
        return {
            expr.name: make_unit_stencil()
        }

    # --- IndexedVar ---
    if isinstance(expr, IndexedVar):
        return {
            expr.name: make_unit_stencil()
        }

    # --- Сложение ---
    if isinstance(expr, Add):
        s1 = to_stencil(expr.left, constants)
        s2 = to_stencil(expr.right, constants)
        return merge_systems(s1, s2)

    # --- Умножение ---
    if isinstance(expr, Mul):
        left, right = expr.left, expr.right
    
        # --- случай: Var * (производная / stencil) ---
        if isinstance(left, Var):
            stencil = to_stencil(right, constants)
            return scale_stencil(stencil, left)   # КЛЮЧЕВОЕ
    
        if isinstance(right, Var):
            stencil = to_stencil(left, constants)
            return scale_stencil(stencil, right)  # КЛЮЧЕВОЕ
    
        # --- константы ---
        if isinstance(left, Const):
            stencil = to_stencil(right, constants)
            return scale_stencil(stencil, left)
    
        if isinstance(right, Const):
            stencil = to_stencil(left, constants)
            return scale_stencil(stencil, right)
    
        # fallback
        s1 = to_stencil(left, constants)
        s2 = to_stencil(right, constants)
        return multiply_systems(s1, s2)
      
    # --- Производная ---
    if isinstance(expr, Derivative):
        inner = expr.expr
        #print(inner)
    
        # 1. базовая переменная
        if isinstance(inner, Var):
            deriv = first_derivative_stencil(expr.axis, expr.scheme)
            s = StencilExpr()
            s.add_term((0,0,0), Const(1))
            return {inner.name: s * deriv}
    
        # 2. уже Derivative → свёртка stencil
        elif isinstance(inner, Derivative):
            # рекурсивно строим stencil для внутренней производной
            inner_stencil = to_stencil(inner, constants)
            outer_deriv = first_derivative_stencil(expr.axis, expr.scheme)
            # свёртка для каждого var
            result = {}
            for var, stencil in inner_stencil.items():
                result[var] = stencil * outer_deriv
            return result
    
        # 3. сложное выражение (Add/Mul)
        else:
            # раскрываем правило дифференцирования через diff
            expanded = diff(inner, expr.axis, expr.scheme, constants)
            return to_stencil(expanded, constants)
          
    if isinstance(expr, TimeDerivative):
        inner = expr.expr
        if isinstance(inner, Var):
            s = StencilExpr()
            s.add_term((0,0,0), Const("1/tau"))
            return {inner.name: s}
    
        raise NotImplementedError("TimeDerivative only for Var")

    raise NotImplementedError(f"Unsupported node in to_stencil: {type(expr)}")

# Константа
def make_const_stencil(c):
    s = StencilExpr()
    s.add_term((0,0,0), c)
    return s
# Единичный stencil (u[i])
def make_unit_stencil():
    s = StencilExpr()
    s.add_term((0,0,0), Const(1))
    return s
# Сложение систем
def merge_systems(s1, s2):
    result = {}
    keys = set(s1) | set(s2)
    for k in keys:
        if k in s1 and k in s2:
            result[k] = s1[k] + s2[k]
        elif k in s1:
            result[k] = s1[k]
        else:
            result[k] = s2[k]

    return result
  
def multiply_systems(s1, s2):
  result = {}

  # Случай: константа * stencil
  for name_l, stencil_l in s1.items():
    for name_r, stencil_r in s2.items():
      # оба __const__
      #print(name_r, name_l)
      if name_l == "__const__" and name_r == "__const__":
          # сохраняем символы + числа
          conv = StencilExpr()
          for offset_l, coef_l in stencil_l.terms.items():
              for offset_r, coef_r in stencil_r.terms.items():
                  offset = tuple(offset_l[i]+offset_r[i] for i in range(3))
                  conv.add_term(offset, Mul(coef_l, coef_r))
          result["__const__"] = conv
          continue

      # const * var
      if name_l == "__const__":
          conv = StencilExpr()
          for offset_l, coef_l in stencil_l.terms.items():
              for offset_r, coef_r in stencil_r.terms.items():
                  offset = tuple(offset_l[i]+offset_r[i] for i in range(3))
                  conv.add_term(offset, Mul(coef_l, coef_r))
          result[name_r] = conv
          continue

      if name_r == "__const__":
          conv = StencilExpr()
          for offset_l, coef_l in stencil_l.terms.items():
              for offset_r, coef_r in stencil_r.terms.items():
                  offset = tuple(offset_l[i]+offset_r[i] for i in range(3))
                  conv.add_term(offset, Mul(coef_l, coef_r))
          result[name_l] = conv
          continue

      # var * var (пока не поддерживается)
      raise NotImplementedError(f"{name_r}*{name_l} пока не поддерживается")

  return result
  
def is_scalar(stencil: StencilExpr):
    return (
        len(stencil.terms) == 1 and
        (0,0,0) in stencil.terms
    )

def get_scalar(stencil: StencilExpr):
    return stencil.terms[(0,0,0)]

# Базовые stencil производных
def first_derivative_stencil(axis, scheme):
  h = f"h_{axis}"
  s = StencilExpr()
  # Центральная разность
  if scheme == Scheme.CENTRAL:
      s.add_term(shift(axis, +1), Const(f"1/(2*{h})"))
      s.add_term(shift(axis, -1), Const(f"-1/(2*{h})"))
      return s
  # Разность вперёд
  if scheme == Scheme.FORWARD:
      s.add_term(shift(axis, +1), Const(f"1/{h}"))
      s.add_term((0,0,0), Const(f"-1/{h}"))
      return s
  # Разность назад
  if scheme == Scheme.BACKWARD:
      s.add_term((0,0,0), Const(f"1/{h}"))
      s.add_term(shift(axis, -1), Const(f"-1/{h}"))
      return s
  raise NotImplementedError("Поддерживаются только центральные разности и разности вперёд-назад")
  
def shift(axis, val):
    if axis == "X":
        return (val,0,0)
    if axis == "Y":
        return (0,val,0)
    if axis == "Z":
        return (0,0,val)
      
def scale_stencil(stencil_dict, coef):
    result = {}
    #print(stencil_dict.keys())
    for var, stencil in stencil_dict.items():
        # создаём копию, чтобы не портить исходный объект
        new_stencil = stencil.copy()

        # масштабируем коэффициенты внутри stencil
        for shift, value in new_stencil.terms.items():
            #print(coef,value)
            new_stencil.terms[shift] = Mul(coef, value)

        result[var] = new_stencil

    return result