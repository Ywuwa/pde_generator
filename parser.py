import ast
from ast_nodes import Var, Const, Add, Mul, Func, TimeDerivative
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

        if func_name.startswith("D") and len(func_name) == 3:
            axis1 = func_name[1]
            axis2 = func_name[2]
            return diff(diff(arg_expr, axis1), axis2)

        if func_name in ["DX", "DY", "DZ"]:
            return diff(arg_expr, func_name[1])
        
        if func_name == "Dt":
            return TimeDerivative(arg_expr)

        if func_name == "L":
            return laplace(arg_expr)

        return Func(func_name, arg_expr)

    raise NotImplementedError("Unsupported syntax")