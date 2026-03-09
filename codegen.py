# codegen.py
import re
from ast_nodes import Expr
from parser import parse_expr
from discretizer import discretize, discretize_ast,generate_cpp
from symbolic import distributive_expand, split_by_variable
from ast_nodes import Add, Mul, IndexedVar
from collections import defaultdict

def decode_offset(index_expr: str):
    """
    index+offset_X-offset_Y -> (1,-1,0)
    """
    dx = dy = dz = 0

    if "offset_X" in index_expr and "+" in index_expr:
        dx = 1
    if "offset_X" in index_expr and "-" in index_expr:
        dx = -1

    if "offset_Y" in index_expr and "+" in index_expr:
        dy = 1
    if "offset_Y" in index_expr and "-" in index_expr:
        dy = -1

    if "offset_Z" in index_expr and "+" in index_expr:
        dz = 1
    if "offset_Z" in index_expr and "-" in index_expr:
        dz = -1

    return dx, dy, dz
  
def offset_to_code(dx, dy, dz):
    def sym(v):
        if v == 1:
            return "R"
        if v == -1:
            return "L"
        return "0"

    return f"{sym(dx)}{sym(dy)}{sym(dz)}"

def flatten_sum(expr):
    if isinstance(expr, Add):
        return flatten_sum(expr.left) + flatten_sum(expr.right)
    return [expr]

def find_indexed_var(expr, var):
    """
    Возвращает (node, index_expr) если найден var[index]
    """
    if isinstance(expr, IndexedVar) and expr.name == var:
        return expr, expr.index
    return None


def extract_coef_and_index(expr, var):

    # случай p[index]
    if isinstance(expr, IndexedVar) and expr.name == var:
        return None, expr.index

    # случай coef * p[index]
    if isinstance(expr, Mul):

        if isinstance(expr.left, IndexedVar) and expr.left.name == var:
            return expr.right, expr.left.index

        if isinstance(expr.right, IndexedVar) and expr.right.name == var:
            return expr.left, expr.right.index

    return None

def collect_stencil_coefficients(lhs_ast, var):

    terms = flatten_sum(lhs_ast)
    coef_map = defaultdict(list)
    for t in terms:
        res = extract_coef_and_index(t, var)
        if not res:
            continue
        coef_ast, index_expr = res
        print(index_expr)
        dx,dy,dz = decode_offset(index_expr)
        coef_map[(dx,dy,dz)].append(coef_ast)

    return coef_map

def generate_stencil_system(lhs_ast: Expr, rhs_ast: Expr, var: str, constants: set):
    coef_map = collect_stencil_coefficients(lhs_ast, var)
    lines = []
    for (dx,dy,dz), coefs in sorted(coef_map.items()):
        print(dx,dy,dz)
        code = offset_to_code(dx,dy,dz)
        coef_expr = " + ".join(generate_cpp(discretize_ast(c, constants)) for c in coefs)
        lines.append(
            f"""const double indexVal_{code} = {coef_expr};
                triplets.emplace_back(index,index+({dx})*offset_X+({dy})*offset_Y+({dz})*offset_Z,indexVal_{code});"""
        )

    lines.append("")
    lines.append(
        f"B[index] = -({generate_cpp(discretize_ast(rhs_ast, constants))});"
    )
    print(lines)
    return lines

def handle_constraint_equation(lhs_ast: Expr, rhs_ast: Expr, var: str, constants_names: set):
    """
    Обрабатывает уравнение вида:
    expression = 0, p

    lhs_ast - левая часть после переноса членов
    rhs_ast - правая часть
    var - переменная (например p)
    """

    # 1. Дискретизация
    lhs_ast = discretize_ast(lhs_ast, constants_names)
    rhs_ast = discretize_ast(rhs_ast, constants_names)

    # 2. Раскрытие скобок
    lhs_ast = distributive_expand(lhs_ast)

    # 3. Генерация stencil
    lines = generate_stencil_system(lhs_ast, rhs_ast, var, constants_names)

    joined_lines = "\n    ".join(lines)
    return joined_lines

def transform_constraint_equation(line: str, constants: set):

    match = re.match(r'(.+)=\s*0\s*,\s*([A-Za-z_]\w*)', line)

    if not match:
        raise ValueError("Invalid constraint equation")

    expr_str = match.group(1).strip()
    variable = match.group(2).strip()

    expr = parse_expr(expr_str)
    lhs, rhs = split_by_variable(expr, variable)
    return handle_constraint_equation(lhs, rhs, variable, constants)



#========================== time equations
def collect_variables(expr, varset, constants: set):
    """
    Рекурсивно собирает переменные в подвыражениях в множество varset
    """
    from ast_nodes import Var, Add, Mul, Func, Derivative

    if isinstance(expr, Var):
        if expr.name not in constants:
            varset.add(expr.name)

    elif isinstance(expr, (Add, Mul)):
        collect_variables(expr.left, varset, constants)
        collect_variables(expr.right, varset, constants)

    elif isinstance(expr, Func):
        collect_variables(expr.arg, varset, constants)

    elif isinstance(expr, Derivative):
        collect_variables(expr.expr, varset, constants)

def handle_time_equation(line: str, constant_names: set):
    """
    Обрабатывает выражение вида:
        Dt(u) + RHS = 0
    Собирает переменные в множество varset.
    Генерирует:
        u1[index] = u[index] + tau*( - RHS )
    """
    match = re.match(r'Dt\(([^)]+)\)\s*\+\s*(.+)\s*=\s*0', line)
    if not match:
        raise ValueError("Неверный формат временного уравнения")

    var = match.group(1).strip()
    rhs_str = match.group(2).strip()

    rhs_expr = parse_expr(rhs_str)

    varset = set()
    collect_variables(rhs_expr, varset, constant_names)
    varset.add(var)

    rhs_cpp = discretize(rhs_expr, constant_names)

    update_code = f"{var}1[index] = {var}[index] + tau*({rhs_cpp});"

    return update_code, sorted(varset)



#=========================== final code generation
def generate_src_n_header(constants: list, equations: list, implicit_eq: list):
    """
    Собирает воедино cpp, hpp файлы
    """
    constant_names = {n for n, _ in constants}
    args = []
    total_code = []
    for equation in equations:
      code, variables = handle_time_equation(
          equation,
          constant_names
      )
      #print(code)
      for v in variables:
          args.append(f"std::vector<double>& {v}")
          args.append(f"std::vector<double>& {v}1")
      
      total_code.append(code)
      
    for eq in implicit_eq:
      code = transform_constraint_equation(eq, constants)
      #print(code)
      total_code.append(code)
      
    args.extend([
        "const uint offset_X",
        "const uint offset_Y",
        "const uint offset_Z",
        "const double h_X",
        "const double h_Y",
        "const double h_Z",
        "const double tau",
        "const size_t dimSize"
    ])

    signature = ",\n               ".join(args)

    head_lines = [f"    const double {n} = {v};" for n, v in constants]
    head_code = "\n    ".join(head_lines)
    equations_code = "\n    ".join(total_code)

    return f"""
#include "../headers/generated.hpp"
void generated({signature})
{{
  {head_code}

  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {{
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {{
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {{
      uint index (k*offset_Z + j*offset_Y + offset_X);
      {equations_code}
      }}
    }}
  }}
}}
""", f"""
#include "settings.hpp"
void generated({signature});
"""