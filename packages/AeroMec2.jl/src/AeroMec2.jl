#
# AeroMec2.jl --
#
# AeroMec.jl is a module that aims at simulating the evolution of some aerosol systems
# under the assumption of the box model (all microphysical quantites are spatially uniform)
#
# All implemented evolution-equations are dimensionless
#
# This module depends on the package Distribution.jl
#
#------------------------------------------------------------------------------
#
# This file is part of the AeroMec module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2021,  Matthew Ozon.
#
#------------------------------------------------------------------------------


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
#

"""
       The GDE in this package are dimensionless and discretized. One iterator computes one time step evolution using the Euler scheme and the Jacobian of that same step can be computed at will.

    The mechanisms included in the model are:
    - coagulation (loss and gain)
    - condensation
    - wall losses
    - nucleation as a boundary condition

    Note that it is possible to add extra source/sink terms quite easily.

    The size discretization of the GDE is computed by integrating the equation over the intervals defined by the discretization of the size space and assuming that the size density and the parameters can be approximated as piecewise constant over the discretization intervals. The assumptions, of course, introduce errors which can be quantified; the quantification of the discretization error is easily tractable in the case of a system that does not undergo coagulation.
    Note: the error due to discretization are fairly small, in practice, compare with the error due to the uncertainty in the growth rate (the convection term of the GDE), especially in the case of fine discretization. This can be estimated by dedicated functions implemented in the package.
"""
module AeroMec2

using Statistics
using LinearAlgebra

# export the functions that implement the coagulation mechanism and its jacobian (w.r.t. the number concentration)
export coagulation!, coagulation_loss!, jacobian_coagulation!, jacobian_coagulation_loss!, coagulation_coefficient!, init_coagulation_loop_indices, init_coagulation_loop_indices! # this last one generates index and gain arrays for the full coagulation so that it runs faster

# export the functions for the condensation
export CondensationGrowth!, jacobian_condensation!

# the losses by wall deposition #LATER: all sorts of linear diagonal damping terms
export WallDeposition!, jacobian_wall_deposition!, WallDepositionRate

# and the nucleation #LATER: all sorts of gain and wheel terms. #NOTE: this is a constant term w.r.t. number concentration, therefore, there is no jacobian!
export Nucleation!

# export the structure that contains the physical parameter of the system
export AeroSys

# export the function that the end user will use
export iter!, jacobian_GDE

export cond_err, loss_err, discretization_err

# define the type AeroSys which contains the mechanisms involved in the system and some characteristic values
include("type.jl")

# coagulation mechanism
include("coagulation.jl")

# condensation mechanism
include("condensation.jl")

# nucleaion mechanism #TODO should be change for source and wheel terms (afine part of the evolution equation)
include("nucleation.jl")

# wall loss mechanism
include("loss.jl")

# iterations: Euler discretisation of the ODE #TODO add the RK4 and FEM? (FEM implemented in AeroSys3.jl)
include("iteration.jl")

end # module
