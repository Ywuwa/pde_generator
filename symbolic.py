from ast_nodes import Const, Var, Add, Mul, Func, Derivative, AXES, Scheme
# calculus
def is_const(expr, constants):
    return isinstance(expr, Var) and expr.name in constants
def is_number(x):
    return isinstance(x, (int, float))
  
def diff(expr, axis, scheme=Scheme.CENTRAL, constants=None):
    """
    Дифференцирует AST-выражение по правилам
    """
    if isinstance(expr, Const):
      return Const(0)
    
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
    
    if isinstance(expr, Mul):  # NEW
        a, b = expr.left, expr.right
        da = diff(a, axis, scheme, constants)
        db = diff(b, axis, scheme, constants)

        # если одна сторона константа → упрощаем
        if constants:
            if is_const(a, constants):
                return Mul(a, db)
            if is_const(b, constants):
                return Mul(da, b)

        # общее правило Лейбница
        return Add(
            Mul(da, b),
            Mul(a, db)
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

    raise NotImplementedError(f"Unsupported differentiation case: {expr}")


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
  """
  Нормализация AST:
  - убирает 0
  - сворачивает константы
  - нормализует Mul
  """

  if isinstance(expr, (Const, Var)):
      return expr

  if isinstance(expr, Add):
      left = simplify_ast(expr.left)
      right = simplify_ast(expr.right)

      # 0 + x
      if isinstance(left, Const) and left.value == 0:
          return right
      if isinstance(right, Const) and right.value == 0:
          return left

      # Const + Const
      if isinstance(left, Const) and isinstance(right, Const):
        if is_number(left.value) and is_number(right.value):
            return Const(left.value + right.value)
        else:
            return Add(left, right)

      return Add(left, right)

  if isinstance(expr, Mul):
      left = simplify_ast(expr.left)
      right = simplify_ast(expr.right)

      # 0 * x
      if (isinstance(left, Const) and left.value == 0) or \
         (isinstance(right, Const) and right.value == 0):
          return Const(0)

      # 1 * x
      if isinstance(left, Const) and left.value == 1:
          return right
      if isinstance(right, Const) and right.value == 1:
          return left

      # Const * Const
      if isinstance(left, Const) and isinstance(right, Const):
        if is_number(left.value) and is_number(right.value):
            return Const(left.value * right.value)
        else:
            # это символическое произведение (например Q1*Q2)
            return Mul(left, right)

      # нормализация: Const всегда слева
      if isinstance(right, Const) and not isinstance(left, Const):
          left, right = right, left

      return Mul(left, right)

  return expr