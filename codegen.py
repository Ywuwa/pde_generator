# codegen.py
import re
from parser import parse_expr
from discretizer import discretize_ast, generate_cpp
from symbolic import distributive_expand, split_by_variable
from ast_nodes import Expr, Const, Add, Mul, IndexedVar
from filegen import generate_equations_hpp, generate_equations_cpp, generate_velocity_residual_cpp, generate_compute_flow_cpp

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

def flatten_sum(expr: Expr) -> [Expr]:
    stack = [expr]
    result = []
    while stack:
        e = stack.pop()
        if isinstance(e, Add):
            stack.append(e.left)
            stack.append(e.right)
        else:
            result.append(e)

    return result

def extract_coef_and_index(expr: Expr, var: str):
    # p[index]
    if isinstance(expr, IndexedVar) and expr.name == var:
        return None, expr.index

    # произведение
    if isinstance(expr, Mul):
        indexed = None
        others = []
        stack = [expr]
        while stack:
            node = stack.pop()
            if isinstance(node, Mul):
                stack.append(node.left)
                stack.append(node.right)
            elif isinstance(node, IndexedVar) and node.name == var:
                indexed = node
            else:
                others.append(node)

        if indexed is None:
            return None

        # собираем коэффициент
        if not others:
            coef = None
        else:
            coef = others[0]
            for o in others[1:]:
                coef = Mul(coef, o)
        
        return coef, indexed.index

    return None

def collect_stencil_coefficients(lhs_ast: Expr, var: str, constants: set):
    terms = flatten_sum(lhs_ast)
    coef_map = defaultdict(list)
    extra_rhs = ""
    for t in terms:
        res = extract_coef_and_index(t, var)
        if not res:
            continue
        coef_ast, index_expr = res
        dx,dy,dz = decode_offset(index_expr)
        check_coef_str = generate_cpp(discretize_ast(coef_ast, constants))
        if (dx,dy,dz)==(0,0,0) and "1/tau" in check_coef_str:
          extra_rhs = check_coef_str + f"*{var}[{index_expr}]"
        coef_map[(dx,dy,dz)].append(coef_ast)

    return coef_map, extra_rhs

def generate_stencil_system(lhs_ast: Expr, rhs_ast: Expr, 
                            var: str, constants: set, ind: int):
    coef_map, extra_rhs = collect_stencil_coefficients(lhs_ast, var, constants)
    extra_rhs_flag = False
    if len(extra_rhs) > 1: extra_rhs_flag = True
    lines = []
    for (dx,dy,dz), coefs in sorted(coef_map.items()):
        code = offset_to_code(dx,dy,dz)
        coef_expr = " + ".join(generate_cpp(discretize_ast(c, constants)) for c in coefs)
        lines.append(
            f"""const double indexVal_{code} = {coef_expr};
                triplets{ind}.emplace_back(index,index+({dx})*offset_X+({dy})*offset_Y+({dz})*offset_Z,indexVal_{code});"""
        )          

    lines.append("")
    lines.append(
        f"B{ind}[index] = {extra_rhs} -({generate_cpp(discretize_ast(rhs_ast, constants))});"
    )
    return lines, extra_rhs_flag

def handle_constraint_equation(lhs_ast: Expr, rhs_ast: Expr, 
                               var: str, constants_names: set, ind: int):
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
    lines, extra_rhs_flag = generate_stencil_system(lhs_ast, rhs_ast, var, constants_names, ind)

    joined_lines = "\n    ".join(lines)
    return joined_lines, extra_rhs_flag

def transform_constraint_equation(line: str, constants: set, index: int):

    match = re.match(r'(.+)=\s*0\s*,\s*([A-Za-z_]\w*)', line)
    if not match:
        raise ValueError("Invalid constraint equation")

    expr_str = match.group(1).strip()
    variable = match.group(2).strip()

    expr = parse_expr(expr_str)
    varset = set()
    collect_variables(expr, varset, constants)
    lhs, rhs = split_by_variable(expr, variable)
    code, extra_rhs_flag = handle_constraint_equation(lhs, rhs, variable, constants, index)
    return code, sorted(varset), extra_rhs_flag, variable



#========================== time equations
def collect_variables(expr, varset: set, constants: set):
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
    # уравнение вида Dt(...) + RHS = 0
    match = re.match(r'Dt\(([^)]+)\)\s*\+\s*(.+)\s*=\s*0', line)
    if match:
      var = match.group(1).strip()
      rhs_str = match.group(2).strip()
  
      rhs_expr = parse_expr(rhs_str)
  
      varset = set()
      collect_variables(rhs_expr, varset, constant_names)
      varset.add(var)
  
      rhs_cpp = generate_cpp(discretize_ast(rhs_expr, constant_names))
  
      update_code = f"{var}1[index] = {var}[index] + tau*({rhs_cpp});"
      residual_code = f"""
      resTerm = {var}1[index] - ( {var}[index] + tau*({rhs_cpp}) );
      vectorResidual += resTerm*resTerm;
      """
  
      return update_code, sorted(varset), residual_code
    # уравнение вида Dt(...) = 0
    match = re.match(r'Dt\(([^)]+)\)\s*=\s*0', line)
    if match:
      var = match.group(1).strip()  
      varset = set()
      collect_variables(Const(0), varset, constant_names)
      varset.add(var)  
      update_code = f"{var}1[index] = {var}[index];"
      residual_code = f"""
      resTerm = {var}1[index] - ( {var}[index] );
      vectorResidual += resTerm*resTerm;
      """ 
      return update_code, sorted(varset), residual_code
    
    raise ValueError("Invalid time equation format") 
      



#=========================== final code generation
def generate_src_n_header(constants: list, equations: list, implicit_eq: list):
    """
    Собирает воедино cpp, hpp файлы
    """
    constant_names = {n for n, _ in constants}
    print(constant_names)
    signature_args = []       # common functions' signature
    time_eq_input_args = []   # arguements of the explicit time equations
    res_input_args = []       # arguements of the velocity residual function
    
    time_eq_code = []         # explicit time equations code
    residual_code = []        # velocity residual code
    
    for equation in equations:
      code, variables, vel_res_code = handle_time_equation(
          equation,
          constant_names
      )
      #print(code)
      for v in variables:
          signature_args.append(f"std::vector<double>& {v}")
          signature_args.append(f"std::vector<double>& {v}1")
          time_eq_input_args.append(f"{v}")
          time_eq_input_args.append(f"{v}1")
          res_input_args.append(f"{v}Exac")
          res_input_args.append(f"{v}")
      
      time_eq_code.append(code)
      residual_code.append(vel_res_code)
      
    impl_signature_args = []
    impl_eq_input_args = []
    impl_eq_code = []
    index = 0
    for eq in implicit_eq:
      code, variables, extra_rhs_flag, var = transform_constraint_equation(eq, constant_names, index)
      for v in variables:
        if v == var and extra_rhs_flag:
          impl_signature_args.append(f"std::vector<double>& {v}")
          impl_eq_input_args.append(f"{v}0")
        elif v == var and extra_rhs_flag == False:
          continue # skip this variable as it is hidden in stencil coefs
        else:
          impl_signature_args.append(f"std::vector<double>& {v}")
          impl_eq_input_args.append(f"{v}")
      impl_signature_args.append(f"std::vector<Eigen::Triplet<double>>& triplets{index}")
      impl_signature_args.append(f"Eigen::VectorXd& B{index}")
      impl_eq_input_args.append(f"triplets{index}")
      impl_eq_input_args.append(f"B{index}")
      impl_eq_code.append(code)
      index += 1
      
    common_signature_extension = [
        "const uint offset_X",
        "const uint offset_Y",
        "const uint offset_Z",
        "const double h_X",
        "const double h_Y",
        "const double h_Z",
        "const double tau",
        "const size_t dimSize"
    ]
    signature_args.extend(common_signature_extension)
    impl_signature_args.extend(common_signature_extension)
    
    common_input_extension = ["offsetX", "offsetY", "offsetZ", 
                                "hX", "hY", "hZ", "tau", "dimSize"]
    time_eq_input_args.extend(common_input_extension)
    impl_eq_input_args.extend(common_input_extension)
    res_input_args.extend(common_input_extension)

    signature = ",\n               ".join(signature_args)
    time_eq_input = ", ".join(time_eq_input_args)    
    residual_input = ", ".join(res_input_args)    
    impl_signature = ",\n               ".join(impl_signature_args)
    impl_eq_input = ", ".join(impl_eq_input_args) 
    

    head_lines = [f"    const double {n} = {v};" for n, v in constants]
    head_code = "\n".join(head_lines)
    time_equations_code = "\n    ".join(time_eq_code)
    impl_equations_code = "\n    ".join(impl_eq_code)
    velocity_residual_code = "\n    ".join(residual_code)
    velocity_residual_cpp = generate_velocity_residual_cpp(signature, head_code, 
                                                           velocity_residual_code)
    compute_flow_cpp_code = generate_compute_flow_cpp(
      time_eq_input, impl_eq_input, residual_input, velocity_residual_cpp)
    generated_cpp = generate_equations_cpp(
      signature, head_code, time_equations_code, impl_signature, impl_equations_code)
    generated_hpp = generate_equations_hpp(signature, impl_signature)

    return generated_cpp, generated_hpp, compute_flow_cpp_code