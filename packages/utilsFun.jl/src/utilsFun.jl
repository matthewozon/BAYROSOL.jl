#
# utilsFun.jl --
#
# utilsFun.jl is the Julia implementation of a few functions that may be often used
#
#------------------------------------------------------------------------------
#
# This file is part of the utilsFun module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2019-2020,  Matthew Ozon.
#
#------------------------------------------------------------------------------

module utilsFun

using Printf

# export the accesible functions: logistic function and (1/a)*log(1+e^{ax})
export softMax,      softMaxA,      softMaxInv,      softMaxInvA
export softMaxDeriv, softMaxDerivA, softMaxInvDeriv, softMaxInvDerivA
export logistic, logisticInv, logisticDeriv
export cauchy, cauchy_deriv

# basic numerical integration
export riemann

# linear basis functions (it can always be useful)
export e_k, e_0, e_M



# a few useful functions
include("utils.jl")
include("riemann.jl")
include("basisFun.jl")



end # module
