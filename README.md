[![AeroMeas CI](https://github.com/matthewozon/BAYROSOL.jl/actions/workflows/CI_AeroMeas.yml/badge.svg)](https://github.com/matthewozon/BAYROSOL.jl/actions/workflows/CI_AeroMeas.yml)

# Retrieval of process rate parameters in the General Dynamic Equation (GDE) for aerosols using Bayesian state estimation

## Packages specifically developed

### EKF
The [EKF.jl](packages/EKF.jl) package offers a framework for the Kalman filter and Fixed Interval Kalman Smoother in their linear and extended forms. It is used for estimating the state of observed time varying system along with their covariances/uncertainty ranges.

### Aerosol mechanisms: GDE
A homogeneous-aerosol-system time evolution may be described by the GDE for aerosols. The package [AeroMec2.jl](packages/AeroMec2.jl) implements a time-and-size discrete solver of the GDE in terms of concentrations; it supports the following mechanisms:
- growth/shrinkage by condensation/evaporation (the advection term)
- nucleation or formation of new small particles at the smaller end of the size discretization spectrum (boundary condition)
- linear losses, e.g. wall or dilution losses (dampening term that is proportional to the size density or concentration)
- coagulation loss and gain (a quadratic term under an integral or sum)
and it allows for an estimation of the errors due to the discretization as well as an estimation of the uncertainty due to eronous values of the parameters, provided that one can estimate a correct upper boundary of the errors in the parameters' values.

### Aerosol measurement: mobility analyzer
An aerosol system is not directly observable, one must observe the system through the eyes of a measurement device to obtain any data. The true state of the system can be retrieved from the data if the measurement model is known (it is a necessary condition, but it might not be sufficient). The [AeroMeas.jl](packages/AeroMeas.jl) package provides some implementations of model measurement device and a more realistic modelization of an SMPS.

### Other packages
Some less specific packages have been developed to make the implementation and the readability easier.
- [utilsFun.jl](packages/utilsFun.jl): implements some often used functions, such as the logistic function.
- [myPlot.jl](packages/myPlot.jl): for creating contour plots which are often used to depict particle size distribution (or other quantities), the piece of code can become a bit cumbersome with native matplotlib functions. The myPlot.jl package makes it a bit easier to plot such quantities by implementing more abstract plotting functions, e.g. imshowData or displayLogData2D.
- [StochProc.jl](packages/StochProc.jl): the idea motivating the implementation of this package was to be able to generate **Stoch**astic **Proc**esses, form a very abstract point of view, by using the order of the process, its dimension and few other parameters. It became a bit more oriented towards covariance computation, more specifically, the computation of the covariance of time varying stochastic processes, e.g. first order vector processes. It is used in the parameter estimation codes, for instance, to encode the size correlation of some parameters.


## Synthetic data generation
The codes for the data simulation is structured into two distinct stages, 1) the simulation of the "true" aerosol system ---solving the GDE for a given set of parameters, e.g. nucleation and growth rates--- and 2) the simulation of the data acquisition ---whose noisiness is controlled by the detector-sample-flow. Two systems are being simulated:
- [NPF event](data_simulation/nucleation_event): in a chamber, some particles reside and are being lost by deposition and coagulation while at some point, for a limited amount of time, some vapor is introduced and induces both New Particle Formation and growth. The measurements of concentration ranges from about 14 nm to 700 nm
- [Steady State](data_simulation/steady_state): in a chamber, some vapor is continuously being injected at a constant rate which triggers nucleation and growth by condensation. The system is observed at the lower end of the particle sizes, almost from nucleation size.

## Acquired data (simulated or from lab experiments)
The purpose of the package being the inversion of the density and the estimation of the GDE's parameters, it became evident that the method should tested against ``real'' data that have been acquired during a lab experiment. The data used to test the method are available in the folders [1952_02](data/1952_02.zip), [1802_01](data/1802_01.zip) and [1906_03](data/1906_03.zip),  and some more simulated data can be found in the folders [0000_00](data/0000_00.zip), [0000_01](data/0000_01.zip), [0000_02](data/0000_02.zip), [0000_03](data/0000_03.zip), [0000_04](data/0000_04.zip).

## Parameter estimation
From the data, either generated by the aforementioned simulations or measured from a lab experiment, it is possible to find the values of the parameters of the GDE that best describe the aerosol system (conditionally to the evolution and measurement model). The two cases [NPF event](parameter_estimation/NE_estimation) and [Steady State](parameter_estimation/SS_estimation) are processed using the Fixed Interval Kalman Smoother, the GDE for aerosols and the measurement models.
And more cases based on DMA-train measurements are treated in [DMA train estimation](parameter_estimation/DMA_train_estimation) also using the FIKS and the GDE.
In all cases, the estimation is run by calling the scripts "main.jl" in the REPL. Please make sure that you have led all the dependencies, set the right LOAD_PATH and change the variables related to the data (e.g. input_folder wich is the path where the data are located).

## DMA-train model
The DMA-train models used for the estimations are given in the folders where the data are.



## General note

Make sure to have the right path and file names before running these codes. The codes are mostly "plug-and-play", especially the packages, but the data simulation and parameter estimation may need some adjustment due to file path differences. The codes may be run with newer versions of Julia, and packages, but there is no guaranty of results.

For the packages installation, nothing more than having the proper dependencies listed in [DEP.md](DEP.md) is needed. To make life easier, and so that the "using" commands works, you need to add the path of each package in the LOAD_PATH variable. To do so, you need to edit the startup.jl file (which might be located at "~/.julia/config/startup.jl") and add the command:

> push!(LOAD_PATH,"/path/to/pkg/pkg.jl/src/")

with the proper path to the package source folder (i.e. the src/ folder of the developed packages). In earlier version of Julia, the startup file may not exists, but instead, the juliarc.jl file does; it has a similar purpose. One might find the file using the command locate startup.jl which list all the instance of files and folder containing the string "startup.jl".

For Windows installations, it can happen that the config folder does not exist, so what you can do is to create a config folder in the julia folder (where julia is installed) and in the config folder, create a startup.jl file and populate it with the appropriate command that needs to be executed when starting up Julia.
For the LOAD_PATH variable you should add in the file:

> push!(LOAD_PATH,"C:\\\\path\\\\to\\\\package\\\\name1.jl\\\\src\\\\")
>
> push!(LOAD_PATH,"C:\\\\path\\\\to\\\\package\\\\name2.jl\\\\src\\\\")

using backslashes instead of forward slashes, where you should modify the string "C:\\\\path\\\\to\\\\package\\\\name[12].jl\\\\src\\\\" so that it points to the location of the sources of the local package.


At some point a proper installation guide will appear, but for now, you can follow the instruction above, or for each package, you may use the add command from Pkg. This way, the packages are directly handled by Julia, no need to update the LOAD_PATH. For instance, for the EKF package
```
import Pkg
Pkg.add(url="https://github.com/matthewozon/BAYROSOL.jl",subdir="packages/EKF.jl")
```


## Refs
For the proof of concept, you may look into [1] where the focus is set on developing the methodology and apply it to simulated data only. An application of the methodology to DMA-train data can be found in [2] where we show that the method does a fairly nice job at estimating -- time varying -- nucleation rates and condensation rates. The method has also been applied in the framework of merging data from different sources, e.g. SMPS DMA and CPC, and compared to other methods in [3]. The starting point of this repository is the zenodo repository [BAYROSOL](https://zenodo.org/record/4450492#.YrGUpjVBzuo)


- [1] M. Ozon et al., Retrieval of process rate parameters in the general dynamic equation for aerosols using Bayesian state estimation: BAYROSOL1.0, GMD, 2021, Vol. 14, p. 3715--3739, [DOI: 10.5194/gmd-14-3715-2021](https://www.doi.org/10.5194/gmd-14-3715-2021)
- [2] M. Ozon et al., Aerosol formation and growth rates from chamber experiments using Kalman smoothing, ACP, 2021, Vol. 21, p. 12595–12611, [DOI: 10.5194/acp-21-12595-2021](https://www.doi.org/10.5194/acp-21-12595-2021)
- [3] D. Stolzenburg et al., Combining instrument inversions for sub-10 nm aerosol number size-distribution measurements, JAS, Vol. 159, p. 105862, [DOI: 10.1016/j.jaerosci.2021.105862](https://doi.org/10.1016/j.jaerosci.2021.105862)
- [4] M. Ozon, D. Stolzenburg and L. Dada, BAYROSOL1.1 [DOI: 10.5281/zenodo.4450492](https://doi.org/10.5281/zenodo.4450492)


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4450492.svg)](https://doi.org/10.5281/zenodo.4450492)
