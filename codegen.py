# codegen.py
import re
from parser import parse_expr
from discretizer import discretize


def collect_variables(expr, varset, constants):
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


def handle_time_equation(line, constant_names):
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


def generate_src_n_header(constants, equations):
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
      print(code)
      for v in variables:
          args.append(f"std::vector<double>& {v}")
          args.append(f"std::vector<double>& {v}1")
      
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