# main.py
import re
from codegen import process_implicit, process_explicit, generate_signature_n_input
from filegen import generate_equations_hpp, generate_equations_cpp, \
                    generate_velocity_residual_cpp, generate_compute_flow_cpp

def parse_input_file(filename):
  constants = []
  explicit_eq = []
  implicit_eq = []

  with open(filename, "r") as f:
    for line in f:
      line = line.strip()
      if not line:
          continue
        
      match_impl = re.match(r'(.+)=\s*0\s*,\s*([A-Za-z_]\w*)', line)
      if match_impl:
          implicit_eq.append((match_impl.group(1), match_impl.group(2)))
          continue
          
      match_expl = re.match(r'Dt\(([^)]+)\)\s*\+\s*(.+)\s*=\s*0', line)
      if match_expl:
          explicit_eq.append((match_expl.group(1), match_expl.group(2)))
          continue

      match = re.match(r'^([A-Za-z_]\w*)\s*=\s*(.+)$', line)
      if match:
          constants.append((match.group(1), match.group(2)))
          continue

      raise ValueError(f"Unknown line: {line}")

  if len(explicit_eq) == 0 and len(implicit_eq) == 0 :
      raise ValueError("No equation found")

  return constants, explicit_eq, implicit_eq
    
if __name__ == "__main__":
    constants, explicit_eq, implicit_eq = parse_input_file("input.txt")
    print(constants)
    print(explicit_eq)
    print(implicit_eq)
    
    cppConst, impl_eq_code = process_implicit(implicit_eq, constants)
    print("GEN:", cppConst + '\n' + impl_eq_code)
    print('\n')
    cppConst, expl_eq_code, velocity_residual_code = \
      process_explicit(explicit_eq, constants)
    print("GEN:", cppConst + '\n' + expl_eq_code)
    
    impl_signature, impl_eq_input, expl_signature, expl_eq_input, residual_input = \
      generate_signature_n_input(implicit_eq, explicit_eq, constants)
      
    generated_cpp = generate_equations_cpp(
      expl_signature, cppConst, expl_eq_code, impl_signature, impl_eq_code)
    
    generated_hpp = generate_equations_hpp(expl_signature, impl_signature)
    
    velocity_residual_cpp = generate_velocity_residual_cpp(
      expl_signature, cppConst, velocity_residual_code)
    
    compute_flow_cpp_code = generate_compute_flow_cpp(
      expl_eq_input, impl_eq_input, residual_input, velocity_residual_cpp)

    with open("src/generated.cpp", "w") as out:
        out.write(generated_cpp)
    with open("headers/generated.hpp", "w") as out:
        out.write(generated_hpp)
    with open("src/compute_flow.cpp", "w") as out:
        out.write(compute_flow_cpp_code)

    print("Done.")