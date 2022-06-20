#
# StochProc.jl --
#
# StochProc.jl is the Julia implementation of some Stochastic Processes
#
#------------------------------------------------------------------------------
#
# This file is part of the StochProc module which is licensed under the MIT "Expat" License:
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

"""
    The idea motivating the implementation of this package was to be able to generate **Stoch**astic **Proc**esses, form a very abstract point of view, by using the order of the process, its dimension and few other parameters. It became a bit more oriented towards covariance computation, more specifically, the computation of the covariance of time varying stochastic processes, e.g. first order vector processes. It is used in the parameter estimation codes, for instance, to encode the size correlation of some parameters.
"""
module StochProc

using LinearAlgebra
using Polynomials #This package may need some modification to make it compatible with the used version of Julia
using StatsBase
using Statistics
using DSP
using ToeplitzMatrices # replaces SpecialMatrices

# because we need those apparently (probably no longer)
import Base: getindex, setindex!, push!, haskey, keys # , start, done, next


# export the stochastic process generators
export  SP_1T, SP_1TC, SP_1TCV, SP_2T, SP_2TC, SP_2TCV, SP_2TCV_mat, SP_1T_V, SP_1TC_V

# export the estimators
export mean_and_varince_estimation, percentile_estimation, histo, multivariate_probability, spherical_coordinate, uniform_sample_unit_sphere, uniform_sample_unit_ball, region_high_probability, percentile_estimation_ni

# covariance computation
export space_covariance_chol, space_covariance_array, space_covariance_range_limit, covariance_process_2nd, covariance_2nd_order_space_1st_order_time, covariance_process_array, covariance_process_range_limit
export covariance_second_order_process

# some very generic way to create random processes
export linSP, iterator, timeBlockMatrix, initGaussRep

# Cholesky factorization wrapped
include("cholWrap.jl")

# covariance generation: size dependence covariance matrix
include("SP_mod.jl")

# test: a unified framework
include("type.jl")

# and the iterator
include("StochProcLin.jl") #TODO: this could be used in the following files in order to have uniform notations

# first order process generation: scalar valued time series
include("StochProc1stOrderScalar.jl")

# first order process generation: scalar valued time series
include("StochProc2ndOrderScalar.jl")

# first order process generation: vector time series
include("StochProc1stOrderVector.jl")

# numerical estimation
include("NE.jl")

end # module

