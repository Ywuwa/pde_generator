# ast_nodes.py
AXES = ["X", "Y", "Z"]
class Scheme:
    CENTRAL = "C"
    FORWARD = "F"
    BACKWARD = "B"
    
class Expr:
    """
    Базовый класс всех узлов выражения.
    """
    pass

class Var(Expr):
    """
    Переменная
    """
    def __init__(self, name):
        self.name = name
    def __str__(self):
      return self.name

class IndexedVar(Expr):
    """
    Индексированная переменная
    """
    def __init__(self,name,index):
        self.name=name
        self.index=index
    def __str__(self):
      return f"{self.name}[{self.index}]"

class Const(Expr):
    """
    Числовая константа
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
      return str(self.value)


class Add(Expr):
    """
    Сложение
    """
    def __init__(self, left, right):
        self.left = left
        self.right = right
    def __str__(self):
      return f"({self.left} + {self.right})"


class Mul(Expr):
    """
    Умножение
    """
    def __init__(self, left, right):
        self.left = left
        self.right = right
    def __str__(self):
      return f"({self.left} * {self.right})"


class Func(Expr):
    """
    Нелинейная функция: sin(u), exp(u), log(u) и т.д.
    """
    def __init__(self, name, arg):
        self.name = name
        self.arg = arg


class Derivative(Expr):
    """
    Производная с указанием схемы:
    C - central
    F - forward
    B - backward
    """
    def __init__(self, expr, axis, scheme=Scheme.CENTRAL):
        self.expr = expr
        self.axis = axis
        self.scheme = scheme
    def __str__(self):
      return f"({self.axis} * {self.scheme})"
        
class TimeDerivative(Expr):
    """
    Производная по времени: Dt(u)
    """
    def __init__(self, expr):
        self.expr = expr