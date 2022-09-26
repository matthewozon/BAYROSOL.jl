# General Dynamic Equation for aerosol

The GDE in this package are dimensionless and discretized. One iterator computes one time step evolution using the Euler scheme and the Jacobian of that same step can be computed at will.

The mechanisms included in the model are:
- coagulation (loss and gain)
- condensation
- wall losses
- nucleation as a boundary condition

Note that it is possible to add extra source/sink terms quite easily.

The size discretization of the GDE is computed by integrating the equation over the intervals defined by the discretization of the size space and assuming that the size density and the parameters can be approximated as piecewise constant over the discretization intervals. The assumptions, of course, introduce errors which can be quantified; the quantification of the discretization error is easily tractable in the case of a system that does not undergo coagulation.
Note: the error due to discretization are fairly small, in practice, compare with the error due to the uncertainty in the growth rate (the convection term of the GDE), especially in the case of fine discretization. This can be estimated by dedicated functions implemented in the package.


# Important note

No significant improvement of this package is expected anymore at this time.


# Dependencies and end-user functions
## Dependencies
### Install

The package dependencies are listed in the file [Project.toml](https://github.com/matthewozon/BAYROSOL.jl/blob/master/packages/AeroMec2.jl/Project.toml).
The installation of the package can be done through the package manager (`import Pkg`) as follows:

`Pkg.add(url="https://github.com/matthewozon/BAYROSOL.jl.git",rev="master",subdir="packages/AeroMec2.jl")`

### Alternative
If the install failed, you can download the package [AeroMec2](https://github.com/matthewozon/BAYROSOL.jl/edit/master/packages/AeroMec2.jl) and install the dependencies.
The module depends on two other packages that can be installed by running the the following commands in the REPL: "import Pkg; Pkg.add("toto")" where "toto" should be the name of the package to install. The packages are:
- [Statistics](https://docs.julialang.org/en/v1/stdlib/Statistics/): is part of the Julia standard library and implements some usual statistical functions, e.g. var.
- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/): is also part of the standard library of Julia; it implements most linear algebra function and wraps function from the optimized BLAS and LAPACK libraries.



## Important functions
Many functions are implemented in this package but only a few are aimed at the end-user who should call only a small subset of function in order to simulate the evolution of an aerosol system. Here is a list of the functions that should be called:

- **AeroSys**: creates an object/structure AeroSys which contains the information required to solve the GDE. It should be called with at least one argument, the array of particle diameters (in meter) arrange in an increasing manner and log-scaled (i.e. d<sub>i</sub> = d<sub>0</sub> r<sup>i-1</sup>) or linear-scaled (i.e. d<sub>i</sub> = d<sub>i-1</sub> + &Delta;<sub>d</sub>). Once a structure is created, one may modify manually the value of some variables (of the structure) such as:
  - beta0: characteristic value of the coagulation coefficients (key: beta_c)
  - GR0: characteristic value of the growth rate (key: GR_c)
  - gamma0: characteristic value of the loss rate (key: gamma_c)
  - J0: characteristic value of the nucleation rate (key: J_c)
  - t0: characteristic time (key: t_c)
  - x0: characteristic concentration (key: x_c)
  - is_coa and is_coa_gain: Boolean variables both set to false by default. If is_coa is set to true, then the loss by coagulation is accounted for in the system, but not the gain. If one wants to simulate the gain by coagulation, one must set both is_coa and is_coa_gain to true; setting is_coa_gain to true and is_coa to false will result in a system without any coagulation mechanism. **Note**: turning on the coagulation mechanism requires an explicit computation of the coefficients (**beta**, a square array) either using the function **coagulation_coefficient!** defined in the file *coagulation.jl* or by computing your own coefficients. **Note**: turning on the coagulation gain means an considerable increase in computation time and memory, as well as an extra explicit initialization of two index arrays and two gain arrays by calling the function **init_coagulation_loop_indices!**.
  - is_con: Boolean variable, set to false by default. If the simulated system must involve growth by condensation, the variable must be set to true
  - is_nuc: Boolean variable, set to false by default. Set to true if nucleation is happening
  - is_los: Set to true is some linear loss mechanism, e.g wall loss or dilution, is accounted for. Note: by linear loss I mean any loss mechanism that can be model as proportional to the concentration (or size density).
An example of creation and initialization of this object is given at the beginning of the file *test/main.jl*.
- **iter!**: the **!** in the syntax mean that the function may modify the arguments it is passed. It computes one time iteration of the aerosol system knowing the current state of the system (by definition, the state of the simulated system is the collection of concentrations, nucleation rate, coagulation coefficients, growth rates and loss rates). The time iteration is solved using a Euler discretization scheme.
- **jacobian_GDE**: this function computes the Jacobian matrix **&delta;** of the discrete time evolution equation of the concentration **N<sub>i</sub><sup>k+1</sup>** = f(**N<sub>i</sub><sup>k</sup>**;J<sup>k</sup>,&beta;<sup>k</sup>,GR<sup>k</sup>,&gamma;<sup>k</sup>), i.e. the matrix whose entries are defined by: <img src="https://latex.codecogs.com/svg.latex?&space;\delta_{i,j}=\frac{\partial N_i^{k+1}}{\partial N_j^{k}}" title="Jacobian element" />
- **discretization_err**: compute the error due to the discretization scheme. It requires an estimation of the first derivative of the density and the knowledge of the growth rate.
- **cond_err**: estimates the error at the current step due to the error in the values of the nucleation rate and the growth rate; the errors in the parameters must be feed to the function.
- **loss_err**: estimates the error in the concentrations due to the error in the loss rates.

All the simulated mechanism are implemented in the files *src/coagulation.jl*, *src/condensation.jl*, *src/loss.jl* and *src/nucleation.jl*.
An example of simulation is given in the file *test/main.jl*; I would suggest to first start trying this simulation without the coagulation mechanism, for the sake of simplicity and clarity.


