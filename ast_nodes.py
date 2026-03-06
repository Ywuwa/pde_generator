# ast_nodes.py
AXES = ["X", "Y", "Z"]

class Expr:
    """
    Базовый класс всех узлов выражения.
    """
    pass


class Var(Expr):
    """
    Переменная (u, v, w,...)
    """
    def __init__(self, name):
        self.name = name


class Const(Expr):
    """
    Числовая константа
    """
    def __init__(self, value):
        self.value = value


class Add(Expr):
    """
    Сложение
    """
    def __init__(self, left, right):
        self.left = left
        self.right = right


class Mul(Expr):
    """
    Умножение
    """
    def __init__(self, left, right):
        self.left = left
        self.right = right


class Func(Expr):
    """
    Нелинейная функция: sin(u), exp(u), log(u) и т.д.
    """
    def __init__(self, name, arg):
        self.name = name
        self.arg = arg


class Derivative(Expr):
    """
    Первая производная произвольного выражения
    """
    def __init__(self, expr, axis):
        self.expr = expr
        self.axis = axis