# main.py
import re
from codegen import generate_src_n_header

def parse_input_file(filename):

    constants = []
    time_equations = []
    implicit_eq = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("Dt("):
                time_equations.append(line)
                continue

            match_eq = re.match(r'(.+)=\s*0\s*,\s*([A-Za-z_]\w*)', line)
            if match_eq:
                implicit_eq.append(line)
                continue

            match = re.match(r'^([A-Za-z_]\w*)\s*=\s*(.+)$', line)
            if match:
                constants.append((match.group(1), match.group(2)))
                continue

            raise ValueError(f"Unknown line: {line}")

    if len(time_equations) == 0 and len(implicit_eq) == 0 :
        raise ValueError("No equation found")

    return constants, time_equations, implicit_eq


if __name__ == "__main__":

    constants, equations, implicit_eq = parse_input_file("input.txt")
    print(equations)
    print(implicit_eq)
    code, header, compute_flow_cpp_code = generate_src_n_header(constants, equations, implicit_eq)

    with open("src/generated.cpp", "w") as out:
        out.write(code)
    with open("headers/generated.hpp", "w") as out:
        out.write(header)
    with open("src/compute_flow.cpp", "w") as out:
        out.write(compute_flow_cpp_code)

    print("Done.")