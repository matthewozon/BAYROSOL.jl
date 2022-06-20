# Measurement operator: Scanning Mobility Particle Sizer

The main goal of this package is to provide a relatively easy-to-create (and use) model of an SMPS. The model is build as a series of module through which the real (input) size density goes and undergoes multiple processes before giving the resulting measurement. Four deterministic and one stochastic processes (or modules) are modeled:
- **impactor**: at the inlet of the SMPS, the big particles are removed by impaction (deterministic)
- **particle charging**: the remaining particles are charged in the charging module. The number of charge per particle depends, among other things, on the size.(deterministic)
- **size selection**: the size classification is in reality operated by a mobility analyzer. (deterministic)
- **detection efficiency**: the size density coming out of the size classifier goes through a counting device whose efficiency may drop drastically for small particles. (deterministic)
- **counting noise**: for optical counting sensors, a Poisson distribution is probably one of the best model. (stochastic)

The impactor and the detection efficiency modules are fairly well modeled by a logistic function which displays a 50% threshold and a sharpness parameter. The measurement noise of the counting device is also well approximated by a Poisson noise with parameter the number of particle passing through the detector. Both the charging of particles and the size selection can not be easily modeled by one function; their model is devised in the mobility space rather than the ill-adapted size space.
The model is based on the user manual of an SMPS (and might not be suited for all SMPS).
Note that the losses are not accounted for --- they could easily be added though.

# Dependencies and end-user functions
## Dependencies
The module depends on several other packages that can be installed by running the the following commands in the REPL: "import Pkg; Pkg.add("toto")" where "toto" should be substitute by the name of the package to install. The packages are:
- [Statistics](https://docs.julialang.org/en/v1/stdlib/Statistics/): is part of the Julia standard library and implements some usual statistical functions, e.g. var.
- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/): is also part of the standard library of Julia; it implements most linear algebra function and wraps function from the optimized BLAS and LAPACK libraries.
- [Distributions](https://juliastats.org/Distributions.jl/stable/starting/): as is said on the official page "The Distributions package provides a large collection of probabilistic distributions and related functions". It is used to generate the counting noise, e.g. Poisson distribution.
- [Printf](https://docs.julialang.org/en/v1/stdlib/Printf/): also part of the standard library, it is used to make string formatting easy.


## Important functions
Many functions are implemented in this package but only a few are aimed at the end-user who should call only a small subset of function in order to simulate the measurement of an aerosol population. Here is a list of the functions that should be called:

- **SMPS3936_transfer_function**: creates a matrix -- the transfer function of the device -- that transforms the size density (or concentrations) into the size density (or concentrations) that reaches the counting sensor. This function must have at least three arguments, the discretization size, the centroid of each channel and a scalar which has no effect at all --- here for compatibility purposes and will probably be dropped in future versions. Additionally, several optional argument ma be passed with the following key words:
  - *s50imp* and *delta50imp*: the impactor 50% threshold and sharpness set by default to 10<sup>-6</sup>m and 0.1x10<sup>-6</sup>m.
  - *s50cpc* and *delta50cpc*: the counting device efficiency 50% threshold and sharpness set by default to 10<sup>-8</sup>m and 10<sup>-9</sup>m.
  - *T* and *Pr*: the temperature and pressure at which the device operates, set by default to 293K and 10<sup>5</sup>Pa.
  - *Nq*: stands for the signed maximal number of chargers acquired by a particle (default value -6)
  - *q_a* and *q_sh* are the aerosol and sheath flux whose default values are 0.3 and 3 L min<sup>-1</sup>
- **SMPS3936**: simulates the measurement given by an SMPS knowing the true size density coming in. Essentially, this function requires the same argument as **SMPS3936_transfer_function** and a few more. To be added at the list of required argument are: the size density of the input population, the discretization lengths (the size difference between two consecutive discretization points), the detector-sample-flow rate and the time spent on one channel (roughly the scanning time devided by the number of channels). Two optional arguments are present but have no effect in this version of the code.
- **DMPS_gaus** and **DMPS_gate**: these function implements fairly poor models of DMPS/SMPS assuming no multiple charging effect, perfect counting efficiency (still with noise) and no impactor. The mobility classifier has a transfer function in the size space which is either gaussian shaped (**DMSP_gaus**) or rectangular (**DMPS_gate**). These function both require the input size distribution, the corresponding discretization size with the interval lengths, the centroids of the channels, the constant ratio between two consecutive channels and the counting volume (volume of sample used for counting in the CPC).

The model is implemented in the files *src/AeroMeas.jl*, *src/charge_probability.jl*, *src/SMPS3936.jl*, *src/SizeDistribution.jl*, *src/DMA.jl* and *src/CPC.jl*.
An example of simulation is given in the file *test/SMPS_model_generator.jl*.

# Bibliography
- Millikan, R., "The general law of fall of a small spherical body through a gas, and its bearing upon the nature of molecular reflection from surfaces", Physical Review, APS, 1923, 22, 1
- Stolzenburg, M., "An ultrafine aerosol size distribution measuring system", University of Minnesota, 1989
- McMurry, P., "The history of condensation nucleus counters", Aerosol Science & Technology, Taylor & Francis, 2000, 33, 297-322
- Flagan, R., "History of electrical aerosol measurements" Aerosol Science and Technology, Taylor & Francis, 1998, 28, 301-380
- Boisdron, Y. & Brock, J., "On the stochastic nature of the acquisition of electrical charge and radioactivity by aerosol particles", Atmospheric Environment (1967), Elsevier, 1970, 4, 35-50
- Wiedensohler, A., "An approximation of the bipolar charge distribution for particles in the submicron size range", Journal of Aerosol Science, Elsevier, 1988, 19, 387-389
- SMPS 3936 TSI user manual.
