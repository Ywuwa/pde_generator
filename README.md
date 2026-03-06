## General info

This is a study pet-project, PDE-generator, which builds C++ code from a Navier-Stocks PDE, described inside a file. The aim of this project is to help simply transfer PDE to code and check their performance on ABC-flow. The domain is a simple cube.

## Python files

+ **main.py** - main file, reads input file, calls the parser, calls generator
+ **codegen.py** - handles equations, collects constants, generates files
+ **parser.py** - parse equations, creates custom AST
+ **discretizer.py** - creates C++ code according to AST
+ **symbolic.py** - describes difference rules
+ **ast_nodes.py** - inherits Python AST-functionality

## C++ files

+ src
  - inout (input-output functionality)
  - main
  - mesh_n_model (init mesh and model parameters, initial conditions, precise solution)
  - compute_flow (construct scheme, solve equations)
  - generated (file, generated via Python PDE_generator)
+ headers
  - settings (contains model parameters structure and includes some basic libraries)
+ config (contains model data aka number of scheme (not needed in current project), domain X-Y-Z size, domain partition, time partition, Reynolds number etc)

## Other
+ input.txt (input file)
+ makefile (builds C++ project)

## HOW_TO_WRITE_PDE
1. Write constants first.
2. If PDE has partial time derivative, write it first (for ex., Dt(u) + ...)
3. If PDE has no partial time derivative, point out which variable is about to be computed (for ex., "... = 0, p" means we had to solve the system "A*p = b")

> Absolute path were set inside settings.hpp, change it in order to be able run program\
> Eigen is needed to run C++ program

