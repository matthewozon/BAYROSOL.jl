#
# EKF.jl --
#
# EKF.jl is the Julia implementation of the Kalman Filter in its Extended version.
# It can perfectly be use for the standard version.
#
#------------------------------------------------------------------------------
#
# This file is part of the EKF module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2018-2019,  Matthew Ozon.
#
#------------------------------------------------------------------------------

module EKF

using LinearAlgebra # contains the pinv function that I use for inverting square matrices (I could also change it to / or inv... but it seems to be better because it is possible to control the inversion of ill-conditioned matrices)

# because we need those apparently
import Base: getindex, setindex!, push!, haskey, keys  #, start, done, next

# the exported functions that need to be implemented (problem specific) (maybe also some constants)
export update_jacobian!, update_model_covariance, update_data_covariance, update_measurement_jacobian!, apply_model!, apply_perfect_measure!, set_jacobian!, set_measurement_jacobian!

# export
export one_step, one_step_fixed_interval, KF, KF_FIS, wsKF

# the type that contains all the variables used in the algorithm
include("type.jl")

# define the function that must be overloaded # TODO: wrapp the functions in a type which is then passed to the Kalman Filters
function update_jacobian! end
function update_measurement_jacobian! end
function update_model_covariance end
function update_data_covariance end

function apply_model! end
function apply_perfect_measure! end
function set_jacobian! end
function set_measurement_jacobian! end

#the function for the temperature estimation. Here, it is required to know the workspace ws,
include("ekf_mod.jl")




end # module
