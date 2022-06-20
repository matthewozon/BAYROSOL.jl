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

#################################################
###      Cholesky factorisation wrapping      ###
#################################################

# compute the Cholesky factor of a symmetric matrix that has negative eigen values because of numerical unstablities or because the eigenvalues are really small in absolute value compare to the largest one.
function chol_PS_L_almost_positive(A::Array{Cdouble,2},tau::Cdouble=1.0e-6)
    if ( (tau<=0.0) | (tau>=1.0) )
        throw("The value of the second argument is not good: tau \$\\in\$ ]0,1[")
    end
    # compute eigen decomposition
    D,P = eigen(Symmetric(A))
    # find the threshold index
    idx_D = sortperm(abs.(D),rev=true)
    Dmax = abs(D[idx_D[1]])
    # idx_star = findfirst(abs(D[idx_D]).<=(tau*Dmax))
    println(D[idx_D])
    println(D[idx_D].<=(tau*Dmax))
    idx_star = findfirst(D[idx_D].<=(tau*Dmax))
    println(idx_star)
    if idx_star==0
        throw("either the eigen values are too negative or the smallest eigenvalue is not small enough")
    end
    if idx_star==1
        throw("The eigenvalues are negative, right?")
    end
    if any(D[idx_D[1:idx_star-1]].<0.0)
        throw("The negative values are to big: the approximation would not be relevant!")
    end
    # replace the tail of the spectrum
    N = length(D)
    if idx_star!=N
        D_star = D[idx_D[idx_star]]
        for i in idx_star:N
            D[idx_D[i]] = D_star
        end
    else
        D[idx_D[idx_star]] = tau*Dmax
    end
    # return
    # full(chol(P*diagm(D)*P',Val{:L}))
    # cholesky(P*diagm(D)*P',Val{true}).L # it keeps track of informations
    cholesky(Symmetric(P*diagm(D)*P')).L
end

# wrap cholesky decomposition for positive semi-definite cases
function chol_PS_L(A::Array{Cdouble,2})
    L = zeros(size(A))
    if any(isnan.(A))
        println("The Cholesky decomposition cannot be compute for a matrix containing NaN. Returning the null matrix")
    else
        if any(isinf.(A))
            println("The Cholesky decomposition cannot be compute for a matrix containing infinite values. Returning the null matrix")
        else
            # try first to compute the Cholesky decomposition for the positive definite case:
            try
                # A = L*L' with L a triangular matrix
                LL = cholesky(Symmetric(A))
                L = convert(Matrix,LL.L)
            catch msgErr
                println(msgErr)
                println("I am now try to compute the Cholesky factor on a positive definite approximation of the matrix: it may not be good!")
                try
                    LLL = chol_PS_L_almost_positive(A)
                    L = convert(Matrix,LLL)
                catch
                    throw(msgErr)
                end
            end
        end
    end
    # return
    L
end
