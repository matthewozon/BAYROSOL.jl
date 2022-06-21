# Julia
## Installation
If you haven't installed Julia on your computer or if you have never used it before, you may want to take a look at the official [Julia website](https://julialang.org).
## Version
Julia Version 1.7.2 (Commit [bf53498635] (2022-02-06 15:21 UTC)) was used for developing the packages and other codes. Here is some extra information for compatibility issues:

- OS: Linux (x86_64-linux-gnu)
- CPU: AMD Ryzen 9 5900HX with Radeon Graphics
- WORD_SIZE: 64
- LIBM: libopenlibm
- LLVM: libLLVM-12.0.1 (ORCJIT, znver3)

Note that [BAYROSOL1.1](https://doi.org/10.5281/zenodo.4450492) was developed with the following configuration:
Julia Version 1.3.1 (Commit [2d5741174c] (2019-12-30 21:36 UTC))

- OS: Linux (x86_64-pc-linux-gnu)
- CPU: Intel(R) Core(TM) i7-7500U CPU @ 2.70GHz
- WORD_SIZE: 64
- LIBM: libopenlibm
- LLVM: libLLVM-6.0.1 (ORCJIT, skylake)

# Development interface
I personally used [Juno](https://junolab.org/) for developing my Julia codes (but not for running them), here is a link to the [installation page](http://docs.junolab.org/latest/man/installation/).
The recommended IDE is now [VS code](https://www.julia-vscode.org/) which I currently use.

# Package dependency
List of the Julia packages upon which the codes depends:

- [PyPlot](https://github.com/JuliaPy/PyPlot.jl): standard plotting library based on matplotlib (version [d330b81b] PyPlot v2.10.0)
- [Printf](https://docs.julialang.org/en/v1/stdlib/Printf/): also part of the standard library, it is used to make string formatting easy.
- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/): is also part of the standard library of Julia; it implements most linear algebra function and wraps function from the optimized BLAS and LAPACK libraries. (BLAS: libopenblas (OpenBLAS 0.3.5  USE64BITINT DYNAMIC_ARCH NO_AFFINITY Haswell MAX_THREADS=16) LAPACK: libopenblas64_)
- [Statistics](https://docs.julialang.org/en/v1/stdlib/Statistics/): is part of the Julia standard library and implements some usual statistical functions, e.g. var.
- [StatsBase](https://juliastats.org/StatsBase.jl/stable/#StatsBase.jl-Documentation-1): from the Julia website "StatsBase.jl is a Julia package that provides basic support for statistics". (version [2913bbd2] StatsBase v0.33.2)
- [Distributions](https://juliastats.org/Distributions.jl/stable/starting/): as is said on the official page "The Distributions package provides a large collection of probabilistic distributions and related functions". It is used to generate the counting noise, e.g. Poisson distribution. ([31c24e10] Distributions v0.24.18)
- [ToeplitzMatrices](https://github.com/JuliaMatrices/ToeplitzMatrices.jl): convenient package for creating special matrices, e.g. Toeplitz. (version [c751599d] ToeplitzMatrices v0.7.0)
- [DataFrames](https://juliadata.github.io/DataFrames.jl/stable/): tabular data manipulation (I use it to load and save data) (version [a93c6f00] DataFrames v1.3.2)
- [Dates](https://docs.julialang.org/en/v1/stdlib/Dates/): date and time manipulation and formatting
- [CSV](https://juliadata.github.io/CSV.jl/stable/index.html): manipulation of csv files (e.g. read, write) (version [336ed68f] CSV v0.10.4)
- [Polynomials](https://juliamath.github.io/Polynomials.jl/stable/): from the website "Julia package that provides basic arithmetic, integration, differentiation, evaluation, and root finding over dense univariate polynomials". This package may need some modification to make it compatible with the used version of Julia. (version [f27b6e38] Polynomials v2.0.17)
- [DSP](https://github.com/JuliaDSP/DSP.jl): is a digital signal processing package which I use for the convolution function. (version [717857b8] DSP v0.7.5)
- [Interpolations](http://juliamath.github.io/Interpolations.jl/latest/): this package implement interpolation and extrapolation so that it can be computed easily. (version [a98d9a8b] Interpolations v0.13.5)

If some of the packages are not installed once Julia has been installed, you may install them manually. If the package is part of the official list of Julia packages, you can run, in the Julia REPL, the following command
> import Pkg; Pkg.add("NameOfThePkg.jl")

If the package is not (yet) an official one, you should follow the instruction given by the developers of the said package. In most case, the installation procedure is described on the website where the package is available.
