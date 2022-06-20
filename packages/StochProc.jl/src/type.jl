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



# a containter that describes the stochastic process

# define the type mSP as the type of possible covariances: scalar or matrix
mSP = Union{Cdouble,Array{Cdouble,2}}
# and the type vSP as the type of possible processes: scalar or vector
vSP = Union{Cdouble,Array{Cdouble,1}}

struct linSP
    # order of the process
    procOrd::Int64
    # number of dimensions
    numDim::Int64

    # time evolution model
    B::mSP
    # covariance of the source
    G_source::mSP
    # lower triangular Cholesky factor of the covariance
    L_source::mSP

    #LATER: from the covariance of the source, it is easy to compute the covariance of the process, it might be interesting to add it to the object creation
    # # covariance of the process
    # G_proc::mSP
    # # lower triangular Cholesky factor of the covariance
    # L_proc::mSP

    # default ctor (it is not really meaningful)
    function linSP() #
        new(0,0,0.0,0.0,0.0)
    end

    # ctor with known order (assume a scalar process)
    function linSP(_ord::Int64)
        if _ord==1
            new(_ord,1,0.0,0.0,0.0)
        else
            new(_ord,1,Array(Cdouble,_ord,_ord),0.0,0.0)
        end
    end

    # ctor with known dimensions
    function linSP(_ord::Int64, _dim::Int64)
        if ((_ord==1) & (_dim==1))
            new(_ord,_dim,0.0,0.0,0.0)
        elseif ((_ord!=1) & (_dim==1))
            new(_ord,_dim,0.0,0.0,0.0)
        else
            new(_ord,_dim,Array(Cdouble,_ord*_dim,_ord*_dim),Array(Cdouble,_dim,_dim),Array(Cdouble,_dim,_dim))
        end
    end

    # ctor: with known dimensions and evolution model
    function linSP(_ord::Int64, _dim::Int64,_B::mSP)
        if isa(_B,Cdouble) # the scalar first order case
            if ((_ord==1) & (_dim==1)) # scalar first order
                new(_ord,_dim,_B,0.0,0.0)
            else
                throw("The evolution model does not match the given dimension and order of the process")
            end
        else
            if size(_B,1)!=size(_B,2)
                throw("The evolution model is not an endomorphism: linSP expects a square matrix")
            end
            if (_ord*_dim)!=size(_B,1)
                throw("The dimension of the evolution matrix does not match: (process order)x(dimension of the process)")
            end
            if ((_ord!=1) & (_dim==1)) # scalar n^{th} order
                new(_ord,_dim,_B,0.0,0.0)
            else  # the actual vectorial case
                new(_ord,_dim,_B,Array(Cdouble,_dim,_dim),Array(Cdouble,_dim,_dim))
            end
        end
    end

    # ctor: with known dimensions, evolution model and source covariance
    function linSP(_ord::Int64, _dim::Int64,_B::mSP,_G_source::mSP)
        if isa(_B,Cdouble) & isa(_G_source,Cdouble) # the scalar first order case
            if ((_ord!=1) | (_dim!=1))
                throw("The evolution model does not match the given dimension and order of the process")
            end
            if _G_source<0.0
                throw("The variance has a negative value")
            end
            new(_ord,_dim,_B,_G_source,sqrt(_G_source))
        elseif isa(_B,Array{Cdouble,2}) & isa(_G_source,Cdouble) # n^{th} order scalar process
            if _G_source<0.0
                throw("The variance has a negative value")
            end
            if size(_B,1)!=size(_B,2)
                throw("The evolution model is not an endomorphism: linSP expects a square matrix")
            end
            if (_ord*_dim)!=size(_B,1)
                throw("The dimension of the evolution matrix does not match: (process order)x(dimension of the process)")
            end
            new(_ord,_dim,_B,_G_source,sqrt(_G_source))
        elseif isa(_B,Array{Cdouble,2}) & isa(_G_source,Array{Cdouble,2}) # the actual vectorial case
            if size(_B,1)!=size(_B,2)
                throw("The evolution model is not an endomorphism: linSP expects a square matrix")
            end
            if size(_G_source,1)!=size(_G_source,2)
                throw("The covariance matrix must be square")
            end
            if (_ord*_dim)!=size(_B,1)
                throw("The dimension of the evolution matrix does not match: (process order)x(dimension of the process)")
            end
            if _dim!=size(_G_source,1)
                throw("The covariance matrix does not match the process dimension")
            end
            L = NaN*Matrix{Cdouble}(I,_dim,_dim)
            try
                L = chol_PS_L(_G_source)
            catch errMsg
                println(errMsg)
                println("The cholesky factor is initialized by NaNs: provided a proper factor before using this objet!")
            end
            new(_ord,_dim,_B,_G_source,L)
        else # first order model for a vectorial process: same root for each element of the vector
            if size(_G_source,1)!=size(_G_source,2)
                throw("The covariance matrix must be square")
            end
            if _dim!=size(_G_source,1)
                throw("The covariance matrix does not match the process dimension")
            end
            L = NaN*Matrix{Cdouble}(I,_dim,_dim)
            try
                L = chol_PS_L(_G_source)
            catch errMsg
                println(errMsg)
                println("The cholesky factor is initialized by NaNs: provided a proper factor before using this objet!")
            end
            new(_ord,_dim,_B*Matrix{Cdouble}(I,_dim,_dim),_G_source,L)
        end
    end

    # cptor
    function linSP(ws::linSP) #
        new(copy(ws.procOrd),copy(ws.numDim),copy(ws.B),copy(ws.G_source),copy(ws.L_source))
    end
end
