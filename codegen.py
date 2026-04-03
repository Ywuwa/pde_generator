# codegen.py
from parser import parse_expr, collect_variables
from stencil_builder import to_stencil

def generate_signature_n_input(
    implicit_eq: list, explicit_eq: list, constants_list_pairs):
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
  
  common_input_extension = ["offsetX", "offsetY", "offsetZ", 
                              "hX", "hY", "hZ", "tau", "dimSize"]
  impl_eq_input_args.extend(common_input_extension)
  expl_eq_input_args.extend(common_input_extension)
  res_input_args.extend(common_input_extension)
  
  impl_signature = ",\n               ".join(impl_signature_args)
  impl_eq_input = ", ".join(impl_eq_input_args) 
  expl_signature = ",\n               ".join(expl_signature_args)
  expl_eq_input = ", ".join(expl_eq_input_args)
  residual_input = ", ".join(res_input_args) 
  
  return impl_signature, impl_eq_input, expl_signature, expl_eq_input, residual_input

def build_constants(constants_list):
    # Заменяет список пар на список строковых значений
    return {name: value for name, value in constants_list}
  
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
      expr = parse_expr(expr_str, constants_list_pairs)
      
      # --- строим stencil ---
      system = to_stencil(expr, constants)
  
      # --- упрощение stencil var ---
      system[var].simplify()
  
      # --- генерация ---
      cppA = generate_cpp(system, var, ind)
      cppB = generate_full_rhs(system, var, ind)
      eqs.append(cppA)
      eqs.append(cppB)
      ind += 1
    
    cppConst = generate_constants(constants)
    impl_eq_code = "\n".join(eqs)
    return cppConst, impl_eq_code

def generate_constants(constants):
    lines = []
    for name, value in constants.items():
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
      expr = parse_expr(expr_str, constants_list_pairs)

      # --- строим stencil ---
      system = to_stencil(expr, constants)
  
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
    
    cppConst = generate_constants(constants)
    expl_eq_code = "\n".join(eqs)
    velocity_residual_code = "\n".join(vel_res_code)
    return cppConst, expl_eq_code, velocity_residual_code
  
def generate_rhs_expl(system):
    """
    Строит правую часть B[index]

    system: результат to_stencil(lhs)

    Логика:
    - знак МИНУС (перенос через =)
    """
    rhs_terms = []
    for name, stencil in system.items():
      expr = stencil_to_cpp_expr(stencil, name)
      rhs_terms.append(f"{expr}")
    if not rhs_terms:
      return "0"
    # перенос в правую часть → знак минус
    return " - (" + " + ".join(rhs_terms) + ")"