## General info

This is a study pet-project, PDE-generator, which builds C++ code from a Navier-Stocks PDE, described inside a file. The aim of this project is to help simply transfer PDE to code and check its performance (residual) on ABC-flow solution. The domain is a simple cube.

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
4. Derivatives: Dt, DX, DY, DZ, DXY, DXZ, etc.

## Simple Examples
+ Constants: `C = 1` -> `const double C = 1`. If there is no equations, the following error will be raised: `ValueError: No equation found`
+ Explicit time equation: `Dt(u) = 0` -> `u1[index] = u[index];`, `Dt(u) + DX(u) = 0` -> `u1[index] = u[index] + tau*((u[index + offset_X] - u[index - offset_X]) / (2*h_X));`. If the format is incorrect, the following error will be raised: `ValueError("Invalid time equation format")` 
+ Implicit equation (Eigen solver is used): `Dt(u) + DX(u) + DZZ(v) = 0, u` ->
    ```
    const double indexVal_L00 = (1/(2*h_X) * -1);
    triplets0.emplace_back(index,index+(-1)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_L00);
    const double indexVal_000 = 1/tau;
    triplets0.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_000);
    const double indexVal_R00 = 1/(2*h_X);
    triplets0.emplace_back(index,index+(1)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_R00);
    
    B0[index] = 1/tau*u[index] -((((v[index + offset_Z + offset_Z] + (-1 * v[index + offset_Z - offset_Z])) + ((-1 * v[index - offset_Z + offset_Z]) + v[index - offset_Z - offset_Z])) * 1/(4*h_Z*h_Z)));
    ```

> Current version of the project fully compiles and works ONLY IF the system has 3 explicit velocity time equations and 1 implicit pressure equation. Otherwise it only generates correct "generated.*" files which can be used the way u like.\
> Absolute path was set inside settings.hpp, change it in order to be able run program\
> Eigen is needed to run C++ program

