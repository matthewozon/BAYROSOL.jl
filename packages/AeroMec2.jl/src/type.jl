#
# type.jl --
#
# type.jl is part of the Module AeroMec and contains the definition of a type that contains
# the relevant information about an aerosol particle system
#
#------------------------------------------------------------------------------
#
# This file is part of the AeroMec module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2021,  Matthew Ozon.
#
#-----------------------------------------------------------------------------

# This file is licensed under the MIT "Expat" License:

# Copyright (c) 2020: Matthew Ozon.

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


"""
    **AeroSys**: creates an object/structure AeroSys which contains the information required to solve the GDE. It should be called with at least one argument, the array of particle diameters (in meter) arrange in an increasing manner and log-scaled (i.e. d<sub>i</sub> = d<sub>0</sub> r<sup>i-1</sup>) or linear-scaled (i.e. d<sub>i</sub> = d<sub>i-1</sub> + &Delta;<sub>d</sub>). Once a structure is created, one may modify manually the value of some variables (of the structure) such as:
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
"""
mutable struct AeroSys

    # bin
    nbin::Int64                    # number of bin
    logS::Bool                     # indicates if the scale is logarithmic (always true for now)
    fixS::Bool                     # indicates if the bin centers can evolve with time (always true for now)
    d::Array{Cdouble,1}            # an array containing the bin diameters (center: mean if logS==false, geometric mean if logS==true) [m]
    d0::Cdouble                    # first bin diameter
    cst_r::Cdouble                 # the diameter constant ratio if logS==true, the diameter constant spacing if logS==false []
    cst_v::Cdouble                 # the volume equivalant of cst_r []

    # coagulation
    beta0::Cdouble                 # a characteristic value of the coagulation coefficients
    beta::Array{Cdouble,2}         # the value of the coagulation coeficients [#^{-1}.cm^3.s^{-1}] #TODO LowerTriangular{Cdouble}

    # condensation
    GR0::Cdouble                   # a characteristic value of the condensational growth rate [m.s^{-1}]
    GR::Cdouble                    # the condensational growth rate for every bins [m.s^{-1}] (we assume that this quantity does not depend on the size of the particle)
    scale_GR:: Array{Cdouble,1}    # the scaling factor for each bin ()

    # nucleation
    J0::Cdouble                    # a characteristic value of the nucleation rate [#.cm^{-3}.s^{-1}]

    # loss rates (possibly due to the walls)
    gamma0::Cdouble                # a characteristic value of the loss rate [s^-1]
    LR::Array{Cdouble,1}           # the loss rate values [s^{-1]]

    # time
    t0::Cdouble                    # a characteristic time of the system [s]

    # concentration
    x0::Cdouble                    # a characteristic number concentration of particles [#m^{-3}]

    # active mechanisms #WARNING: no coagulation implemented for linear scales
    is_coa::Bool                   # true if the coagulation is active #NOTE: coagulation_coefficient!(ws::AeroSys) must be called if true, and you may want to change the units using the scaling factor ws.beta0
    is_con::Bool                   # true if the condensation is active
    is_nuc::Bool                   # true if the nucleation is active
    is_los::Bool                   # true if the linear losses are active
    #WARNING: if true, init_coagulation_loop_indices!(ws::AeroSys) must be called prior to running a simulation and the computation time may explode.
    is_coa_gain::Bool              # true is the coagulation involves both loss and gain, false is only losses

    # making coagulation gain a bit easier
    count_coag_index::Int64        # number of index pairs
    Ic::Array{Int64,2}             # index pairs
    a_gain::Array{Cdouble,1}       # relative gain

    # and the same for the jacobian of the coagulation term
    count_coag_indexj::Int64       # number of index pairs
    Icj::Array{Int64,2}            # index pairs
    a_gainj::Array{Cdouble,1}      # relative gain

    # default ctor (it is not really meaningful)
    function AeroSys() #
        new(0,true,true,Array{Cdouble,1}(undef,0),0.0,0.0,0.0, # dimensions
            0.0, Array{Cdouble,2}(undef,0,0),                  # coagulation
            0.0, 0.0, Array{Cdouble,1}(undef,0),               # condensation
            0.0,                                               # nucleation
            0.0, Array{Cdouble,1}(undef,0),                    # losses
            1.0, 1.0,                                          # general normalization constant
            false,false,false,false,false,                     # mechanisms making the system evolve
            0,Array{Int64,2}(undef,0,0),Array{Cdouble,1}(undef,0),
            0,Array{Int64,2}(undef,0,0),Array{Cdouble,1}(undef,0))
    end

    # ctor with known bin centers
    function AeroSys(_d::Array{Cdouble,1};beta_c::Cdouble=0.0, GR_c::Cdouble=0.0, J_c::Cdouble=0.0, gamma_c::Cdouble=0.0, t_c::Cdouble=1.0, x_c::Cdouble=1.0, scale::String="log")
        nbin_ = length(_d)
        logS_ = (scale=="log")
        if logS_
            cst_r_ = mean(_d[2:end]./_d[1:end-1])
            cst_v_ = cst_r_^3
            scale_GR_ = (cst_r_.^(1.5.-collect(1.0:nbin_)))/(cst_r_-1.0) # d0/delta_i
        else # assume linear scale for now. WARNING: no coagulation implemented yet with linear scale
            cst_r_ = 1.0
            cst_v_ = 1.0
            scale_GR_ = (_d[1]/mean(_d[2:end]-_d[1:end-1]))*ones(Cdouble,nbin) # d0/delta
        end
        new(nbin_,logS_,true,_d,_d[1],cst_r_,cst_v_,   # dimensions
            beta_c, zeros(Cdouble,nbin_,nbin_),       # coagulation
            GR_c, 0.0, scale_GR_,                     # condensation
            J_c,                                      # nucleation
            gamma_c, zeros(Cdouble,nbin_),            # losses
            t_c, x_c,                                 # general normalization constant
            false,false,false,false,false,            # mechanisms making the system evolve
            0,Array{Int64,2}(undef,0,0),Array{Cdouble,1}(undef,0),
            0,Array{Int64,2}(undef,0,0),Array{Cdouble,1}(undef,0))
    end

    # cptor
    function AeroSys(ws::AeroSys) #
        new(copy(ws.nbin), copy(ws.logS), copy(ws.fixS), copy(ws.d), copy(ws.d0), copy(ws.cst_r), copy(ws.cst_v),
            copy(ws.beta0), copy(ws.beta),
            copy(ws.GR0), copy(ws.GR), copy(ws.scale_GR),
            copy(ws.J0),
            copy(ws.gamma0), copy(ws.LR),
            copy(ws.t0), copy(ws.x0),
            copy(ws.is_coa), copy(ws.is_con), copy(ws.is_nuc), copy(ws.is_los),copy(ws.is_coa_gain),
            copy(ws.count_coag_index),  copy(ws.Ic),  copy(ws.a_gain),
            copy(ws.count_coag_indexj), copy(ws.Icj), copy(ws.a_gainj))
    end
end
