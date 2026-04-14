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
3. Click on the Run button or press `Shift + Enter` to run the code.

### In Terminal:
1. Navigate to the folder containing your Julia files using the `cd` command.
2. Start Julia by typing `julia` and pressing `Enter`.
3. Run your script by typing `include("your_script.jl")` (replace `your_script.jl` with your actual script name).
