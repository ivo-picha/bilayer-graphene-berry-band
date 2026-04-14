# Description
Julia code for part **b)** of exercise **2** of **HW06** for the _Quantum Transport & Topological Insulators 2026_ course. Calculate band structure, Berry curvature and Chern number of AB-stacked bilayer graphene with and without inversion symmetry.

# Julia Installation Guide

## Installing Julia
1. Download Julia from the [official website](https://julialang.org/downloads/).
2. Follow the installation instructions for your operating system.

## Installing Required Packages
To use the code in this repository, you need to install the following Julia packages:

1. **Plots**
   ```julia
   using Pkg
   Pkg.add("Plots")
   ```

2. **LinearAlgebra** (this package is included in Julia's standard library, so no need to install)

3. **Measures**
   ```julia
   using Pkg
   Pkg.add("Measures")
   ```

## Running the Code
You can run the Julia scripts from either Visual Studio Code or the terminal:

### In Visual Studio Code:
1. Open the folder containing the Julia files.
2. Open the desired Julia file.
3. Click on the Run button or press `Alt + Enter` to run the code.

### In Terminal:
1. Navigate to the folder containing your Julia files using the `cd` command.
2. Run the code in Julia by typing `julia bilayer_graphene_topo.jl` and pressing `Enter`.
