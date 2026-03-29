import ast
from ast_nodes import Var, Const, Add, Mul, Func, TimeDerivative, Scheme
from symbolic import diff, laplace

def parse_expr(expr_str: str, constants: set):
    """
    Парит выражение средствами AST-библиотеки Python.
    Вызывает построение кастомного AST-дерева
    """
    expr_str = expr_str.strip()
    tree = ast.parse(expr_str, mode='eval').body
    return build_ast(tree, constants)

def build_ast(node, constants):
  """
  Построение кастомного AST-дерева
  """
  if isinstance(node, ast.Name):
      return Var(node.id)

  if isinstance(node, ast.Constant):
      return Const(node.value)

  if isinstance(node, ast.BinOp):

      if isinstance(node.op, ast.Add):
          return Add(build_ast(node.left, constants), build_ast(node.right, constants))

      if isinstance(node.op, ast.Sub):
          return Add(build_ast(node.left, constants),
                     Mul(Const(-1), build_ast(node.right, constants)))

      if isinstance(node.op, ast.Mult):
          return Mul(build_ast(node.left, constants), build_ast(node.right, constants))

  if isinstance(node, ast.Call):
    func_name = node.func.id
    arg_expr = build_ast(node.args[0], constants)

    if func_name == "Dt":
        return TimeDerivative(arg_expr)
    # --- Пространственные производные---
    if func_name.startswith("D"):
        name = func_name[1:].replace("_", "")
    
        ops = parse_derivative_chain(name)
    
        result = arg_expr
        for axis, scheme in ops:
            result = diff(result, axis, scheme)
    
        return result

    if func_name == "L":
        return laplace(arg_expr, scheme=Scheme.CENTRAL, constants=constants)

    return Func(func_name, arg_expr)

  raise NotImplementedError("Unsupported syntax")
    

def parse_derivative_chain(name):
    """
    "XBDYF" -> [("X","B"), ("Y","F")]
    "XY"    -> [("X","C"), ("Y","C")]
    "X"     -> [("X","C")]
    """
    tokens = []
    i = 0
    while i < len(name):
        axis = name[i]
        if axis not in "XYZ":
            raise ValueError(f"Invalid axis in derivative: {name}")

        # указан ли тип схемы после указания оси
        if i + 1 < len(name) and name[i+1] in ["F", "B", "C"]:
            scheme = name[i+1]
            i += 2
        else:
            scheme = "C"
            i += 1
        tokens.append((axis, scheme))
    return tokens