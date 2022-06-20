#
# AeroMeas.jl --
#
# AeroMeas.jl is a module that aims at simulating/modelling measurements on an aerosol system
#
# This module depends on the package Distribution.jl
#
#------------------------------------------------------------------------------
#
# This file is part of the AeroMeas module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2021,  Matthew Ozon.
#
#------------------------------------------------------------------------------


# This file is licensed under the MIT "Expat" License:

# Copyright (c) 2021: Matthew Ozon.

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


module AeroMeas

# import the packages used in this package
using Distributions # for the CPC model which uses Poisson distribution for the counting noise
using Statistics
using LinearAlgebra
using Printf # for the message thrown if an exception arises

# export the function that the end user will use
export DMPS_gaus, DMPS_gate # WARNING: these are very bad models of measurement devices
# these should not be used
export DMA_gaus, DMA_gate, DMA_time_avg_gaus, DMA_time_avg_gate, ch_eff, DMA_gaus_OP, chanel_efficiency, chanel_efficiency_gaus
export CPC, CPC_gaus
# better model (only the two following function should be used, the above ones should be ignored)
export SMPS3936_transfer_function, SMPS3936 # intent to model the SMPS that we have in the lab

# define the basics of particle charging
include("charge_probability.jl")

# implement an SMPS model based on its instruction manual: it's already a quite "sofisticated" model
include("SMPS3936.jl")


# the measurement device: this depends on the package Distribution.jl # This is a very simplified model of an SMPS
include("SizeDistribution.jl")


# implement the inversion method described by Kenneth Wolfenbarger 1990 (J. Aerosol Sci., Vol. 21, No. 2, p. 227-247)
#TODO


end # module
