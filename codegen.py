# codegen.py
<<<<<<< HEAD
from parser import parse_expr, collect_variables
from stencil_builder import to_stencil

def generate_signature_n_input(
    implicit_eq: list, explicit_eq: list, constants_list_pairs: list):
  """
  На основе имеющихся уравнений генерирует сигнатуры и входные значения функций
  """
  index = 0
  varset = set()
  
  # сигнатура и входящие аргументы неявных уравнений
  impl_signature_args = []
  impl_eq_input_args = []
  for eq_str in implicit_eq:
    expr_str, var = eq_str
    expr = parse_expr(expr_str, constants_list_pairs)
    # --- превращаем список пар константа-значение в список констант ---
    constants = build_constants(constants_list_pairs)
    collect_variables(expr, varset, constants)
    variables = sorted(varset)
    for v in variables:
        impl_signature_args.append(f"std::vector<double>& {v}")
        impl_eq_input_args.append(f"{v}")
    impl_signature_args.append(
      f"std::vector<Eigen::Triplet<double>>& triplets{index}")
    impl_signature_args.append(f"Eigen::VectorXd& B{index}")
    impl_eq_input_args.append(f"triplets{index}")
    impl_eq_input_args.append(f"B{index}")
    index += 1
  
  # сигнатура и входящие аргументы явных уравнений  
  varset.clear()
  expl_signature_args = []
  expl_eq_input_args = []
  res_input_args = []
  for eq_str in explicit_eq:  
    var, rhs_str = eq_str
    rhs_expr = parse_expr(rhs_str, constants_list_pairs)
    constants = build_constants(constants_list_pairs)
    collect_variables(rhs_expr, varset, constants)
    varset.add(var)
    variables = sorted(varset)
    for v in variables:
      expl_signature_args.append(f"std::vector<double>& {v}")
      expl_signature_args.append(f"std::vector<double>& {v}1")
      expl_eq_input_args.append(f"{v}")
      expl_eq_input_args.append(f"{v}1")
      res_input_args.append(f"{v}Exac")
      res_input_args.append(f"{v}")
      
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
  impl_signature_args.extend(common_signature_extension)
  expl_signature_args.extend(common_signature_extension)
=======
import re
from parser import parse_expr
from discretizer import discretize_ast, generate_cpp, simplify
from symbolic import distributive_expand, split_by_variable
from ast_nodes import Expr, Const, Add, Mul, Var, IndexedVar, Func, Derivative
from filegen import generate_equations_hpp, generate_equations_cpp, generate_velocity_residual_cpp, generate_compute_flow_cpp

from collections import defaultdict

def decode_offset(index_expr: str):
    """
    index+offset_X-offset_Y -> (1,-1,0)
    """
    dx = dy = dz = 0
    index_expr = index_expr.replace(" ", "")

    pattern = r'([+-])(\d*)\*?offset_X'
    matchv = re.match(pattern, index_expr[5:])
    if matchv:
      sign = matchv.group(1)
      num = matchv.group(2)
      number = num if num else "1"
      dx = int(sign + number)

    pattern = r"([+-])(\d*)\*?offset_Y"
    matchv = re.match(pattern, index_expr[5:])
    if matchv:
      sign = matchv.group(1)
      num = matchv.group(2)
      number = num if num else "1"
      dy = int(sign + number)

    pattern = r"([+-])(\d*)\*?offset_Z"
    matchv = re.match(pattern, index_expr[5:])
    if matchv:
      sign = matchv.group(1)
      num = matchv.group(2)
      number = num if num else "1"
      dz = int(sign + number)

    return dx, dy, dz
  
def offset_to_code(dx, dy, dz):
    def sym(v):
        if v > 0:
            return "R"
        if v < 0:
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
        print((dx,dy,dz))
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
    lhs_ast = simplify(lhs_ast)
    rhs_ast = simplify(rhs_ast)
    lhs_ast = discretize_ast(lhs_ast, constants_names)
    rhs_ast = discretize_ast(rhs_ast, constants_names)
    lhs_ast = simplify(lhs_ast)
    rhs_ast = simplify(rhs_ast)
    
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

    expr = parse_expr(expr_str, constants)
    varset = set()
    collect_variables(expr, varset, constants)
    lhs, rhs = split_by_variable(expr, variable)
    code, extra_rhs_flag = handle_constraint_equation(lhs, rhs, variable, constants, index)
    return code, sorted(varset), extra_rhs_flag, variable

#========================== 
# --- Time equations --
#========================== 

def collect_variables(expr, varset: set, constants: set):
    """
    Рекурсивно собирает переменные в подвыражениях в множество varset
    """
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
  
      rhs_expr = parse_expr(rhs_str, constant_names)
      rhs_expr = simplify(rhs_expr)
      varset = set()
      collect_variables(rhs_expr, varset, constant_names)
>>>>>>> backup_branch
  
  common_input_extension = ["offsetX", "offsetY", "offsetZ", 
                              "hX", "hY", "hZ", "tau", "dimSize"]
  impl_eq_input_args.extend(common_input_extension)
  expl_eq_input_args.extend(common_input_extension)
  res_input_args.extend(common_input_extension)
  
<<<<<<< HEAD
  impl_signature = ",\n               ".join(impl_signature_args)
  impl_eq_input = ", ".join(impl_eq_input_args) 
  expl_signature = ",\n               ".join(expl_signature_args)
  expl_eq_input = ", ".join(expl_eq_input_args)
  residual_input = ", ".join(res_input_args) 
  
  return impl_signature, impl_eq_input, expl_signature, expl_eq_input, residual_input

def build_constants(constants_list):
    # Заменяет список пар на список строковых значений
    return {name for name, value in constants_list}
  
def process_implicit(implicit_eq: str, constants_list_pairs: list):
    """
    Обрабатывает:
    Dt(u/v/w) + EXPR(u,v,w) = 0, u/v/w
    """
    # --- превращаем список пар константа-значение в список констант ---
    constants = build_constants(constants_list_pairs)
    
    # --- проходим по всем уравнениям ---
    eqs = []
    ind = 0
    for eq_str in implicit_eq:
      # --- парсинг строки ---
      expr_str, var = eq_str
      expr = parse_expr(expr_str, constants)
      
      # --- строим stencil ---
      system = to_stencil(expr, constants)
      print(system)
  
      # --- упрощение stencil var ---
      system[var].simplify()
  
      # --- генерация ---
      cppA = generate_cpp(system, var, ind)
      cppB = generate_full_rhs(system, var, ind)
      eqs.append(cppA)
      eqs.append(cppB)
      ind += 1
    
    cppConst = generate_constants(constants_list_pairs)
    impl_eq_code = "\n".join(eqs)
    return cppConst, impl_eq_code

def generate_constants(constants):
    lines = []
    for name, value in constants:
        lines.append(f"const double {name} = {value};")
    return "\n".join(lines)  

def decode(offset):
    m = {-1:"L",0:"0",1:"R"}
    return "".join(m[d] for d in offset)

def generate_cpp(system, var: str, ind: int):
    """
    Генерация triplets (матрица A)
    """
    stencil = system[var]
    lines = []
    for offset, coef in stencil.terms.items():
        name = f"indexVal_{decode(offset)}"
        dx,dy,dz = offset
        lines.append(f"const double {name} = {coef};")
        lines.append(
            f"triplets{ind}.emplace_back(index,"
            f"index+({dx})*offset_X+({dy})*offset_Y+({dz})*offset_Z,"
            f"{name});"
        )

    return "\n".join(lines)
  
def generate_full_rhs(system, var: str, ind: int):
    """
    Полная генерация RHS
    """
    rhs = generate_rhs(system, var)
    if "1/tau" in system[var].terms.items():
      rhs = add_time_term(rhs, var)
    return f"B{ind}[index] = {rhs};"

def generate_rhs(system, var):
    """
    Строит правую часть B[index]

    system: результат to_stencil(lhs)
    var: переменная, относительно которой решаем уравнение

    Логика:
    - все stencil НЕ var → в RHS
    - знак МИНУС (перенос через =)
    """
    rhs_terms = []
    for name, stencil in system.items():
      if name == var:
          continue  # это идёт в матрицу A
      expr = stencil_to_cpp_expr(stencil, name)
      rhs_terms.append(f"{expr}")
    if not rhs_terms:
      return "0"
    # перенос в правую часть → знак минус
    return " - (" + " + ".join(rhs_terms) + ")"

def add_time_term(rhs_expr, var):
    """
    Добавляет вклад Dt(u):
    Dt(u) → 1/tau * u[index] переносится в RHS
    """
    time_part = f"(1/tau * {var}[index])"
    if rhs_expr == "0":
        return time_part

    return f"{rhs_expr} + {time_part}"
  
def stencil_to_cpp_expr(stencil, var_name):
    """
    Преобразует stencil в C++ выражение:
    sum coef * var[index + shift]
    
    Например:
    DX(v) → (v[i+1] - v[i]) / h
    """
    terms = []
    for (dx, dy, dz), coef in stencil.terms.items():
        idx = index_shift(dx, dy, dz)
        if var_name == "__const__":
            # константа
            terms.append(f"{coef}")
        else:
            terms.append(f"({coef} * {var_name}[{idx}])")
    return " + ".join(terms) if terms else "0"

def index_shift(dx, dy, dz):
    """
    Генерирует C++ выражение смещения индекса
    """
    parts = ["index"]
    if dx != 0:
        parts.append(f"{dx}*offset_X")
    if dy != 0:
        parts.append(f"{dy}*offset_Y")
    if dz != 0:
        parts.append(f"{dz}*offset_Z")

    return " + ".join(parts)

def process_explicit(explicit_eq: str, constants_list_pairs: list):
    """
    Обрабатывает:
    Dt(u/v/w) + EXPR(u,v,w) = 0
    """
    # --- превращаем список пар константа-значение в список констант ---
    constants = build_constants(constants_list_pairs)   
    
    # --- проходим по всем уравнениям ---
    eqs = []
    vel_res_code = []
    ind = 0
    for eq_str in explicit_eq:
      # --- парсинг строки ---
      var, expr_str = eq_str
      print(expr_str)
      expr = parse_expr(expr_str, constants)

      # --- строим stencil ---
      system = to_stencil(expr, constants)
      print(system)
  
      # --- упрощение stencil var ---
      system[var].simplify()
  
      # --- генерация ---
      common_gen_rhs = generate_rhs_expl(system)
      cppB = f"{var}1[index] = {var}[index] + tau*({common_gen_rhs});"
      residual_code = f"""
      resTerm = {var}1[index] - ( {var}[index] + tau*({common_gen_rhs}) );
      vectorResidual += resTerm*resTerm;
      """
      eqs.append(cppB)
      vel_res_code.append(residual_code)
      ind += 1
    
    cppConst = generate_constants(constants_list_pairs)
    expl_eq_code = "\n".join(eqs)
    velocity_residual_code = "\n".join(vel_res_code)
    return cppConst, expl_eq_code, velocity_residual_code
  
def generate_rhs_expl(system):
=======
      update_code = f"{var}1[index] = {var}[index] - tau*({rhs_cpp});"
      residual_code = f"""
      resTerm = {var}1[index] - ( {var}[index] - tau*({rhs_cpp}) );
      vectorResidual += resTerm*resTerm;
      """
      return update_code, sorted(varset), residual_code, var
    
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
      
#=============================== 
# --- final code generation ---
#=============================== 

def generate_src_n_header(constants: list, equations: list, implicit_eq: list):
>>>>>>> backup_branch
    """
    Строит правую часть B[index]

    system: результат to_stencil(lhs)

    Логика:
    - знак МИНУС (перенос через =)
    """
<<<<<<< HEAD
    rhs_terms = []
    for name, stencil in system.items():
      expr = stencil_to_cpp_expr(stencil, name)
      rhs_terms.append(f"{expr}")
    if not rhs_terms:
      return "0"
    # перенос в правую часть → знак минус
    return " - (" + " + ".join(rhs_terms) + ")"
=======
    constant_names = {n for n, _ in constants}
    #print(constant_names)
    signature_args = []       # common functions' signature
    time_eq_input_args = []   # arguements of the explicit time equations
    res_input_args = []       # arguements of the velocity residual function
    
    time_eq_code = []         # explicit time equations code
    residual_code = []        # velocity residual code
    total_variables = set()
    
    for equation in equations:
      code, variables, vel_res_code, var = handle_time_equation(
          equation,
          constant_names
      )
      for v in variables:
        if v not in total_variables:
          signature_args.append(f"std::vector<double>& {v}")
          time_eq_input_args.append(f"{v}")
          res_input_args.append(f"{v}Exac")
        total_variables.add(v)
      signature_args.append(f"std::vector<double>& {var}1")
      time_eq_input_args.append(f"{var}1")
      res_input_args.append(f"{var}")
        
      time_eq_code.append(code)
      residual_code.append(vel_res_code)
      
    impl_signature_args = []
    impl_eq_input_args = []
    impl_eq_code = []
    index = 0
    total_variables = set()
    for eq in implicit_eq:
      code, variables, extra_rhs_flag, var = transform_constraint_equation(eq, constant_names, index)
      for v in variables:
        if v not in total_variables:
          if v == var and extra_rhs_flag:
            impl_signature_args.append(f"std::vector<double>& {v}")
            impl_eq_input_args.append(f"{v}0")
          elif v == var and extra_rhs_flag == False:
            continue # skip this variable as it is hidden in stencil coefs
          else:
            impl_signature_args.append(f"std::vector<double>& {v}")
            impl_eq_input_args.append(f"{v}")
        total_variables.add(v)
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
        "const size_t dimSize",
        "const model_data& params"
    ]
    signature_args.extend(common_signature_extension)
    impl_signature_args.extend(common_signature_extension)
    
    common_input_extension = ["offsetX", "offsetY", "offsetZ", 
                                "hX", "hY", "hZ", "tau", "dimSize", "params"]
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
>>>>>>> backup_branch
