import ast
from ast_nodes import Var, Const, Add, Mul, Func, TimeDerivative, Scheme
from symbolic import diff, laplace

def parse_expr(expr_str):
    """
    Парит выражение средствами AST-библиотеки Python.
    Вызывает построение кастомного AST-дерева
    """
    expr_str = expr_str.strip()
    tree = ast.parse(expr_str, mode='eval').body
    return build_ast(tree)


def build_ast(node):

    if isinstance(node, ast.Name):
        return Var(node.id)

    if isinstance(node, ast.Constant):
        return Const(node.value)

    if isinstance(node, ast.BinOp):

        if isinstance(node.op, ast.Add):
            return Add(build_ast(node.left), build_ast(node.right))

        if isinstance(node.op, ast.Sub):
            return Add(build_ast(node.left),
                       Mul(Const(-1), build_ast(node.right)))

        if isinstance(node.op, ast.Mult):
            return Mul(build_ast(node.left), build_ast(node.right))

    if isinstance(node, ast.Call):
      func_name = node.func.id
      arg_expr = build_ast(node.args[0])
  
      # --- Производные ---
      if func_name.startswith("D"):
          name = func_name[1:]  # убираем D
  
          # случай DXY, DXX и т.д.
          if len(name) == 2 and all(c in "XYZ" for c in name):
              axis1, axis2 = name
              return diff(diff(arg_expr, axis1), axis2)
  
          # DX, DY, DZ
          if len(name) == 1 and name in "XYZ":
              return diff(arg_expr, name, scheme=Scheme.CENTRAL)
  
          # DXF, DYB и т.д.
          if len(name) == 2 and name[1] in ["F", "B", "C"]:
              return diff(arg_expr, name[0], scheme=name[1])
  
      if func_name == "Dt":
          return TimeDerivative(arg_expr)
  
      if func_name == "L":
          return laplace(arg_expr)
  
      return Func(func_name, arg_expr)

    raise NotImplementedError("Unsupported syntax")