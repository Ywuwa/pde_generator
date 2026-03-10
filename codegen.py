# codegen.py
import re
from ast_nodes import Expr
from parser import parse_expr
from discretizer import discretize, discretize_ast, generate_cpp
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

def collect_stencil_coefficients(lhs_ast: Expr, var: str):
    terms = flatten_sum(lhs_ast)
    coef_map = defaultdict(list)
    for t in terms:
        res = extract_coef_and_index(t, var)
        if not res:
            continue
        coef_ast, index_expr = res
        dx,dy,dz = decode_offset(index_expr)
        coef_map[(dx,dy,dz)].append(coef_ast)

    return coef_map

def generate_stencil_system(lhs_ast: Expr, rhs_ast: Expr, 
                            var: str, constants: set, ind: int):
    coef_map = collect_stencil_coefficients(lhs_ast, var)
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
        f"B{ind}[index] = -({generate_cpp(discretize_ast(rhs_ast, constants))});"
    )
    return lines

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
    lines = generate_stencil_system(lhs_ast, rhs_ast, var, constants_names, ind)

    joined_lines = "\n    ".join(lines)
    return joined_lines

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
    return handle_constraint_equation(lhs, rhs, variable, constants, index), sorted(varset)



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
    residual_code = f"""
    resTerm = {var}1[index] - ( {var}[index] + tau*({rhs_cpp}) );
    vectorResidual += resTerm*resTerm;
    """

    return update_code, sorted(varset), residual_code



#=========================== final code generation
def generate_src_n_header(constants: list, equations: list, implicit_eq: list):
    """
    Собирает воедино cpp, hpp файлы
    """
    constant_names = {n for n, _ in constants}
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
      code, variables = transform_constraint_equation(eq, constants, index)
      for v in variables:
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
    head_code = "\n    ".join(head_lines)
    time_equations_code = "\n    ".join(time_eq_code)
    impl_equations_code = "\n    ".join(impl_eq_code)
    velocity_residual_code = "\n    ".join(residual_code)
    velocity_residual_cpp = generate_velocity_residual_cpp(signature, head_code, 
                                                           velocity_residual_code)
    compute_flow_cpp_code = generate_compute_flow_cpp(
      time_eq_input, impl_eq_input, residual_input, velocity_residual_cpp)

    return f"""
#include "../headers/generated.hpp"
void generated_time_eq({signature})
{{
  {head_code}

  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {{
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {{
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {{
      uint index (k*offset_Z + j*offset_Y + offset_X);
      {time_equations_code}
      }}
    }}
  }}
}}

void generated_impl_eq({impl_signature})
{{
  {head_code}

  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {{
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {{
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {{
      uint index (k*offset_Z + j*offset_Y + offset_X);
      {impl_equations_code}
      }}
    }}
  }}
}}
""", f"""
#include "settings.hpp"
void generated_time_eq({signature});
void generated_impl_eq({impl_signature});
""", compute_flow_cpp_code


def generate_velocity_residual_cpp(signature, head_code, velocity_residual_code):
  return f"""
/*
<var>   - exact value at knot
<var>1  - estimated solution at knot
*/
double velocity_residual({signature})
{{
  {head_code}
  double vectorResidual (0.0);
  //---------------------------- inner knots --------------------------------
  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {{
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {{
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {{
        uint index (k*offsetZ + j*offsetY + i);
        double resTerm (0.0); // residual term
        //! Insert precise values to the scheme, take the difference with the estimated values
        {velocity_residual_code}
      }}
    }}
  }}
  //-------------------------------------------------------------------------
  return std::sqrt(vectorResidual);
}}
"""
  

def generate_compute_flow_cpp(
    time_eq_input: str, impl_eq_input: str, residual_input: str, velocity_residual_cpp: str):
  return f""" 
#include "../headers/inout.hpp"
#include "../headers/mesh_n_model.hpp"
#include "../headers/compute_flow.hpp"
#include "../headers/generated.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

//======================================= FLOW COMPUTATION ========================================
void compute_cube(
  const model_data& params, 
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p0)
{{
  uint tick = 1; // number of time step
  const auto dimSize ( params.domainPartition );  // 1-dimension size
  const uint offsetX = 1;
  const uint offsetY = dimSize + 1;               // j+1 component (i+1 component is just [index+1])
  const uint offsetZ = (dimSize+1) * (dimSize+1); // k+1 component

  const double tau ( params.duration / params.timePartition ); // time step
  const double hX (params.xLen / dimSize ); // x-step
  const double hY (params.yLen / dimSize ); // y-step
  const double hZ (params.zLen / dimSize ); // z-step

  std::ofstream outputFile(params.PATH + params.PATH_log, std::ios::out);
  std::ofstream outputResidualFile(params.PATH + params.PATH_residual, std::ios::out);
  const std::string outputFuncFile = params.PATH;
  // vector size
  const size_t vecSize = (dimSize + 1) * (dimSize + 1) * (dimSize + 1);
  std::vector<double> u1(vecSize, 0.0);
  std::vector<double> v1(vecSize, 0.0);
  std::vector<double> w1(vecSize, 0.0);
  Eigen::VectorXd p(vecSize);
  for (size_t i = 0; i < vecSize; i++) p[i] = p0[i];

  while (tick < params.timePartition + 1)
  {{
    //! velocity exact
    std::vector<double> uExac(vecSize, 0.0);
    initialConditions(uExac, 0, params, tick*tau);
    std::vector<double> vExac(vecSize, 0.0);
    initialConditions(vExac, 1, params, tick*tau);
    std::vector<double> wExac(vecSize, 0.0);
    initialConditions(wExac, 2, params, tick*tau);
    std::vector<double> pExac(vecSize, 0.0);
    initialConditions(pExac, 3, params, tick*tau);

    //! velocity compute
    //---------------------------- inner knots --------------------------------
    generated_time_eq({time_eq_input});
    //-------------------------------------------------------------------------

    //---------------------------- border knots -------------------------------
    //! XY-plane Z = 0 / Z = MAX; Neiman's condition dw/dz = 0
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {{
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {{
        uint index1 (j*offsetY + i);
        uint index2 (dimSize*offsetZ + j*offsetY + i);
        u1[index1] = uExac[index1];
        u1[index2] = uExac[index2];
        v1[index1] = vExac[index1];
        v1[index2] = vExac[index2];

        w1[index1] = wExac[index1];
        w1[index2] = wExac[index2];
      }}
    }}
    //! XZ-plane Y = 0 / Y = MAX; Neiman's condition dv/dy = 0
    for (size_t k = 1; k < dimSize; k++)    // Z-Axis
    {{
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {{
        uint index1 (k*offsetZ + i);
        uint index2 (k*offsetZ + dimSize*offsetY + i);
        w1[index1] = wExac[index1];
        w1[index2] = wExac[index2];
        u1[index1] = uExac[index1];
        u1[index2] = uExac[index2];

        v1[index1] = vExac[index1];
        v1[index2] = vExac[index2];
      }}
    }}
    //! YZ-plane X = 0 / X = MAX; Neiman's condition du/dx = 0
    for (size_t k = 1; k < dimSize; k++)    // Z-Axis
    {{
      for (size_t j = 1; j < dimSize; j++)  // X-Axis
      {{
        uint index1 (k*offsetZ + j*offsetY);
        uint index2 (k*offsetZ + j*offsetY + dimSize);
        w1[index1] = wExac[index1];
        w1[index2] = wExac[index2];
        v1[index1] = vExac[index1];
        v1[index2] = vExac[index2];

        u1[index1] = uExac[index1];
        u1[index2] = uExac[index2];
      }}
    }}
    //-------------------------------------------------------------------------

    //---------------------------- edge knots ---------------------------------
    //! Y = Z = (0 || MAX); dv/dy = 0; dw/dz = 0
    for (size_t i = 1; i < dimSize; i++)
    {{
      // Y = Z = 0
      uint index1 (i);
      // Y = MAX, Z = 0
      uint index2 (dimSize*offsetY + i);
      // Y = 0, Z = MAX
      uint index3 (dimSize*offsetZ + i);
      // Y = MAX, Z = MAX
      uint index4 (dimSize*offsetZ + dimSize*offsetY + i);

      u1[index1] = uExac[index1];
      u1[index2] = uExac[index2];
      u1[index3] = uExac[index3];
      u1[index4] = uExac[index4];

      v1[index1] = vExac[index1];
      v1[index2] = vExac[index2];
      v1[index3] = vExac[index3];
      v1[index4] = vExac[index4];

      w1[index1] = wExac[index1];
      w1[index2] = wExac[index2];
      w1[index3] = wExac[index3];
      w1[index4] = wExac[index4];
    }}
    //! X = Z = (0 || MAX); du/dx = 0, dw/dz = 0
    for (size_t j = 1; j < dimSize; j++)
    {{
      // X = Z = 0
      uint index1 (j*offsetY);
      // X = MAX, Z = 0
      uint index2 (j*offsetY + dimSize);
      // X = 0, Z = MAX
      uint index3 (dimSize*offsetZ + j*offsetY);
      // X = MAX, Z = MAX
      uint index4 (dimSize*offsetZ + j*offsetY + dimSize);

      v1[index1] = vExac[index1];
      v1[index2] = vExac[index2];
      v1[index3] = vExac[index3];
      v1[index4] = vExac[index4];

      u1[index1] = uExac[index1];
      u1[index2] = uExac[index2];
      u1[index3] = uExac[index3];
      u1[index4] = uExac[index4];

      w1[index1] = wExac[index1];
      w1[index2] = wExac[index2];
      w1[index3] = wExac[index3];
      w1[index4] = wExac[index4];
    }}
    //! X = Y = (0 || MAX)' du/dx = 0, dv/dy = 0
    for (size_t k = 1; k < dimSize; k++)
    {{
      // X = Y = 0
      uint index1 (k*offsetZ);
      // X = MAX, Y = 0
      uint index2 (k*offsetZ + dimSize);
      // X = 0, Y = MAX
      uint index3 (k*offsetZ + dimSize*offsetY);
      // X = MAX, Y = MAX
      uint index4 (k*offsetZ + dimSize*offsetY + dimSize);

      w1[index1] = wExac[index1];
      w1[index2] = wExac[index2];
      w1[index3] = wExac[index3];
      w1[index4] = wExac[index4];

      u1[index1] = uExac[index1];
      u1[index2] = uExac[index2];
      u1[index3] = uExac[index3];
      u1[index4] = uExac[index4];

      v1[index1] = vExac[index1];
      v1[index2] = vExac[index2];
      v1[index3] = vExac[index3];
      v1[index4] = vExac[index4];
    }}
    //-------------------------------------------------------------------------

    //------------------------- cube vertices ---------------------------------
    // X = Y = Z = 0
    uint index (0);
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = MAX, Y = Z = 0
    index = dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Z = 0, Y = MAX
    index = dimSize*offsetY;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Y = 0, Z = MAX
    index = dimSize*offsetZ;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Y = MAX, Z = 0
    index = dimSize*offsetY + dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Z = MAX, Y = 0
    index = dimSize*offsetZ + dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // Y = Z = MAX, X = 0
    index = dimSize*offsetZ + dimSize*offsetY;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Y = Z = MAX
    index = dimSize*offsetZ + dimSize*offsetY + dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];

    //-------------------------------------------------------------------------
    // move data from upper time layer
    u.clear(); u = std::move(u1);
    v.clear(); v = std::move(v1);
    w.clear(); w = std::move(w1);
    u1.resize(vecSize); v1.resize(vecSize); w1.resize(vecSize);
    //-------------------------------------------------------------------------
    
    //-------------------------------------------------------------------------
    // barrier, sync point
    //-------------------------------------------------------------------------

    //! initialization of equation entities (Ax = b)
    Eigen::SparseMatrix<double> A(vecSize, vecSize);
    Eigen::VectorXd B0(vecSize);
    std::vector<Eigen::Triplet<double>> triplets0; // entities for filling a sparse matrix

    //! pressure compute
    //---------------------------- inner knots --------------------------------
    generated_impl_eq({impl_eq_input})
    //-------------------------------------------------------------------------

    //---------------------------- border knots -------------------------------
    //! XY-plane Z = 0,1 / Z = MAX-1,MAX; Neiman's condition dp/dz = 0
    for (size_t j = 2; j < dimSize-1; j++)    // Y-Axis
    {{
      for (size_t i = 2; i < dimSize-1; i++)  // X-Axis
      {{
        uint index1 (j*offsetY + i);
        uint index2 (dimSize*offsetZ + j*offsetY + i);
        uint index3 (offsetZ + j*offsetY + i);
        uint index4 ((dimSize-1)*offsetZ + j*offsetY + i);
        // matrix construct via EIGEN triplets 
        // behind low border
        triplets0.emplace_back(index3,index3, 1.0);
        // low border
        triplets0.emplace_back(index1,index1, 1.0);
        // behind up border
        triplets0.emplace_back(index4,index4, 1.0);
        // up border
        triplets0.emplace_back(index2,index2, 1.0);
        B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];
      }}
    }}
    //! XZ-plane Y = 0,1 / Y = MAX-1,MAX; Neiman's condition dp/dy = 0
    for (size_t k = 2; k < dimSize-1; k++)    // Z-Axis
    {{
      for (size_t i = 2; i < dimSize-1; i++)  // X-Axis
      {{
        uint index1 (k*offsetZ + i);
        uint index2 (k*offsetZ + dimSize*offsetY + i);
        uint index3 (k*offsetZ + offsetY + 1);
        uint index4 (k*offsetZ + (dimSize-1)*offsetY + i);
        // matrix construct via EIGEN triplets 
        // behind Y0 border
        triplets0.emplace_back(index3,index3, 1.0);
        // low border
        triplets0.emplace_back(index1,index1, 1.0);
        // behind up border
        triplets0.emplace_back(index4,index4, 1.0);
        // up border
        triplets0.emplace_back(index2,index2, 1.0);
        B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];
      }}
    }}
    
    //! YZ-plane X = 0,1 / X = MAX-1,MAX; Neiman's condition dp/dx = 0
    for (size_t k = 2; k < dimSize-1; k++)    // Z-Axis
    {{
      for (size_t j = 2; j < dimSize-1; j++)  // X-Axis
      {{
        uint index1 (k*offsetZ + j*offsetY);
        uint index2 (k*offsetZ + j*offsetY + dimSize);
        uint index3 (k*offsetZ + j*offsetY + 1);
        uint index4 (k*offsetZ + j*offsetY + dimSize-1);
        // matrix construct via EIGEN triplets (flow out)
        // behind X0 border
        triplets0.emplace_back(index3,index3, 1.0);
        // low border
        triplets0.emplace_back(index1,index1, 1.0);
        // behind up border
        triplets0.emplace_back(index4,index4, 1.0);
        // up border
        triplets0.emplace_back(index2,index2, 1.0);
        B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];
      }}
    }}
    //-------------------------------------------------------------------------

    //---------------------------- edge knots ---------------------------------
    //! Y = Z = (0,1 || MAX,MAX-1); 
    for (size_t i = 2; i < dimSize-1; i++)
    {{
      // Y = Z = 0
      uint index1 (i);
      triplets0.emplace_back(index1,index1, 1.0);
      // Y = MAX, Z = 0
      uint index2 (dimSize*offsetY + i);
      triplets0.emplace_back(index2,index2, 1.0);
      // Y = 0, Z = MAX
      uint index3 (dimSize*offsetZ + i);
      triplets0.emplace_back(index3,index3, 1.0);
      // Y = MAX, Z = MAX
      uint index4 (dimSize*offsetZ + dimSize*offsetY + i);
      triplets0.emplace_back(index4,index4, 1.0);

      B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];

      // Y = Z = 1
      index1 = offsetZ + offsetY + i;
      triplets0.emplace_back(index1,index1, 1.0);
      // Y = MAX-1, Z = 1
      index2 = (dimSize-1)*offsetY + offsetZ + i;
      triplets0.emplace_back(index2,index2, 1.0);
      // Y = 1, Z = MAX-1
      index3 = (dimSize-1)*offsetZ + offsetY + i;
      triplets0.emplace_back(index3,index3, 1.0);
      // Y = MAX-1, Z = MAX-1
      index4 = (dimSize-1)*offsetZ + (dimSize-1)*offsetY + i;
      triplets0.emplace_back(index4,index4, 1.0);
      
      B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];
    }}
    //! X = Z = (0 || MAX); du/dx = 0, dw/dz = 0
    for (size_t j = 2; j < dimSize-1; j++)
    {{
      // X = Z = 0
      uint index1 (j*offsetY);
      triplets0.emplace_back(index1,index1, 1.0);
      // X = MAX, Z = 0
      uint index2 (j*offsetY + dimSize);
      triplets0.emplace_back(index2,index2, 1.0);
      // X = 0, Z = MAX
      uint index3 (dimSize*offsetZ + j*offsetY);
      triplets0.emplace_back(index3,index3, 1.0);
      // X = MAX, Z = MAX
      uint index4 (dimSize*offsetZ + j*offsetY + dimSize);
      triplets0.emplace_back(index4,index4, 1.0);
      
      B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];

      // X = Z = 1
      index1 = offsetZ + j*offsetY + 1;
      triplets0.emplace_back(index1,index1, 1.0);
      // X = MAX-1, Z = 1
      index2 = offsetZ + j*offsetY + dimSize-1;
      triplets0.emplace_back(index2,index2, 1.0);
      // X = 1, Z = MAX-1
      index3 = (dimSize-1)*offsetZ + j*offsetY + 1;
      triplets0.emplace_back(index3,index3, 1.0);
      // X = MAX-1, Z = MAX-1
      index4 = (dimSize-1)*offsetZ + j*offsetY + dimSize-1;
      triplets0.emplace_back(index4,index4, 1.0);
      
      B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];
    }}
    //! X = Y = (0 || MAX)' du/dx = 0, dv/dy = 0
    for (size_t k = 2; k < dimSize-1; k++)
    {{
      // X = Y = 0
      uint index1 (k*offsetZ);
      triplets0.emplace_back(index1,index1, 1.0);
      // X = MAX, Y = 0
      uint index2 (k*offsetZ + dimSize);
      triplets0.emplace_back(index2,index2, 1.0);
      // X = 0, Y = MAX
      uint index3 (k*offsetZ + dimSize*offsetY);
      triplets0.emplace_back(index3,index3, 1.0);
      // X = MAX, Y = MAX
      uint index4 (k*offsetZ + dimSize*offsetY + dimSize);
      triplets0.emplace_back(index4,index4, 1.0);
      
      B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];

      // X = Y = 1
      index1 = k*offsetZ + offsetY + 1;
      triplets0.emplace_back(index1,index1, 1.0);
      // X = MAX-1, Y = 1
      index2 = k*offsetZ + offsetY + dimSize-1;
      triplets0.emplace_back(index2,index2, 1.0);
      // X = 1, Y = MAX-1
      index3 = k*offsetZ + (dimSize-1)*offsetY + 1;
      triplets0.emplace_back(index3,index3, 1.0);
      // X = MAX-1, Y = MAX-1
      index4 = k*offsetZ + (dimSize-1)*offsetY + dimSize-1;
      triplets0.emplace_back(index4,index4, 1.0);
      
      B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];
    }}
    //-------------------------------------------------------------------------

    //------------------------- cube vertices ---------------------------------
    // X = Y = Z = 1
    index = offsetZ + offsetY + 1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = MAX-1, Y = Z = 1
    index = offsetZ + offsetY + dimSize-1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Z = 1, Y = MAX-1
    index = offsetZ + (dimSize-1)*offsetY + 1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Y = 1, Z = MAX-1
    index = (dimSize-1)*offsetZ + offsetY + 1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Y = MAX-1, Z = 1
    index = offsetZ + (dimSize-1)*offsetY + dimSize-1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Z = MAX-1, Y = 1
    index = (dimSize-1)*offsetZ + offsetY + dimSize-1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // Y = Z = MAX-1, X = 1
    index = (dimSize-1)*offsetZ + (dimSize-1)*offsetY + 1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Y = Z = MAX-1
    index = (dimSize-1)*offsetZ + (dimSize-1)*offsetY + dimSize-1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];

    // X = Y = Z = 0
    index = 0;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = MAX, Y = Z = 0
    index = dimSize;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Z = 0, Y = MAX
    index = dimSize*offsetY;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Y = 0, Z = MAX
    index = dimSize*offsetZ;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Y = MAX, Z = 0
    index = dimSize*offsetY + dimSize;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Z = MAX, Y = 0
    index = dimSize*offsetZ + dimSize;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // Y = Z = MAX, X = 0
    index = dimSize*offsetZ + dimSize*offsetY;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Y = Z = MAX
    index = dimSize*offsetZ + dimSize*offsetY + dimSize;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // barrier, sync point
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // call Eigen solver for Pressure (idk, does is support concurrency or not)
    //-------------------------------------------------------------------------
    A.setFromTriplets (triplets0.begin(), triplets0.end());
    // biconjugate gradient stabilized algorithm
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver(A);
    if (solver.info() != Eigen::Success)
    {{
      outputFile << "Can not build preconditioner" << std::endl;
      std::cout << "Can not build preconditioner" << std::endl;
      return;
    }}
    Eigen::VectorXd pHat(vecSize);
    pHat = solver.solveWithGuess(B0, p);
    if (solver.info() != Eigen::Success)
    {{
      outputFile << "Failed to solve the system with Eigen, tick = " << tick << std::endl;
      std::cout << "Failed to solve the system with Eigen, tick = " << tick << std::endl;
      return;
    }}
    //-------------------------------------------------------------------------
    // refresh Pressure values (transfer Eigen-vector to vector)
    //-------------------------------------------------------------------------

    for (size_t i = 0; i < vecSize; ++i) {{p0[i] = pHat[i]; p[i] = pHat[i];}}

    // residual
    //-------------------------------------------------------------------------
    const double velResidual = velocity_residual({residual_input});
    outputResidualFile << std::scientific << velResidual << std::endl;
    //-------------------------------------------------------------------------

    funcOutput(outputFuncFile, "/v1", std::to_string(tick), ".txt", u, params, false);
    funcOutput(outputFuncFile, "/v2", std::to_string(tick), ".txt", v, params, false);
    funcOutput(outputFuncFile, "/v3", std::to_string(tick), ".txt", w, params, false);
    funcOutput(outputFuncFile, "/p", std::to_string(tick), ".txt", p0, params, false);
    tick += 1;
    if ((tick % 100 == 0)) std::cout << "tick: " << tick << '\n';
  }}
  for (size_t i = 0; i < vecSize; ++i) p0[i] = p[i];
  std::cout << "final tick: " << tick << '\n';
}}
//=================================================================================================



//======================================== RESIDUALS ==============================================
{velocity_residual_cpp}
//=================================================================================================

"""