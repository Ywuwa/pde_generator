# main.py
import re
from codegen import generate_src_n_header

def parse_input_file(filename):

    constants = []
    equations = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("Dt("):
                equations.append(line)
                continue

            match = re.match(r'^([A-Za-z_]\w*)\s*=\s*(.+)$', line)
            if match:
                constants.append((match.group(1), match.group(2)))
                continue

            raise ValueError(f"Unknown line: {line}")

    if len(equations) == 0:
        raise ValueError("No Dt equation found")

    return constants, equations


if __name__ == "__main__":

    constants, equations = parse_input_file("input.txt")
    print(equations)
    code, header = generate_src_n_header(constants, equations)

    with open("src/generated.cpp", "w") as out:
        out.write(code)
    with open("headers/generated.hpp", "w") as out:
        out.write(header)

    print("Done.")