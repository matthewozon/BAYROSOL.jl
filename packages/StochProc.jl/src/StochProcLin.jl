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



function iterator(x::vSP,A::linSP)
    if A.numDim==1 # scalar
        if A.procOrd==1 # first order
            # compute iteration
            z = A.B*x .+ A.L_source*randn()
        else # n^{th} order
            z = A.B*x # deterministic evolution
            z[1] = z[1] .+ A.L_source*randn() # addition of the random source
        end
    else # vector
        if A.procOrd==1 # first order
            z = A.B*x .+ A.L_source*randn(A.numDim)
        else # n^{th} order
            eta = A.L_source*randn(A.numDim)
            z = A.B*x
            z[1:A.procOrd:end] = z[1:A.procOrd:end] .+ eta
        end
    end
    z
end

function timeBlockMatrix(B_sub::mSP,dim::Int64)
    if isa(B_sub,Cdouble)
        B = B_sub*Matrix{Cdouble}(I,dim,dim)
    else
        if size(B_sub,1)!=size(B_sub,2)
            throw("The time evolution model is expected to be an endomorphism")
        end
        ord = size(B_sub,1)
        B = zeros(Cdouble,ord*dim,ord*dim)
        for i in 1:dim
            Rg = (i-1)*ord+1:i*ord
            B[Rg,Rg] = B_sub
        end
    end
    B
end

function initGaussRep(mu::vSP,L::mSP,rep::Int64)
    if (isa(mu,Cdouble) & isa(L,Cdouble))
        X = mu .+ L*randn(rep)
    elseif (isa(mu,Array{Cdouble,1}) & isa(L,Cdouble))
        dim = length(mu)
        X = Array{Cdouble,1}(undef,rep*dim)
        for i in 1:rep
            X[i:rep:end] = mu .+ L*randn(dim)
        end
    elseif (isa(mu,Cdouble) & isa(L,Array{Cdouble,2}))
        if size(L,1)!=size(L,2)
            throw("The noise model is expected to be a square matrix")
        end
        dim = size(L,1)
        X = Array{Cdouble,1}(undef,rep*dim)
        for i in 1:rep
            X[i:rep:end] = mu .+ L*randn(dim)
        end
    else
        if size(L,1)!=size(L,2)
            throw("The noise model is expected to be a square matrix")
        end
        if size(L,1)!=length(mu)
            throw("The mean of the noise has not the same dimmension as the covariance")
        end
        dim = length(mu)
        X = Array{Cdouble,1}(undef,rep*dim)
        for i in 1:rep
            X[i:rep:end] = mu .+ L*randn(dim)
        end
    end
    X
end
