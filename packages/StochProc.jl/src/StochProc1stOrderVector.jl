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



###############################################################################
###   vectorial  processes                                                  ###
###############################################################################


# non-linear first order vectorial process
function SP_1T_V(f_time::Function,L_source::Array{Cdouble,2},Nt::Int64, x0::Array{Cdouble,1}) # L_source could be made LowerTriangular or Diagonal or Cdouble
    Nd = length(x0)
    # allocate an array
    X = Array{Cdouble,2}(undef,Nd,Nt)
    # initi the time series
    X[:,1] = x0
    # generate the time series
    for i in 2:Nt
        X[:,i] = f_time(X[:,i-1]) + L_source*randn(Nd)
    end
    # return
    X
end

# linear first order vectorial process: using a scalar as the evolution model (same evolution model for every elements)
function SP_1T_V(r_time::Cdouble,L_source::Array{Cdouble,2},Nt::Int64, x0::Array{Cdouble,1}) # L_source could be made LowerTriangular or Diagonal or Cdouble
    function f_time(X::Array{Cdouble,1}) r_time*X end
    SP_1T_V(f_time,L_source,Nt, x0)
end

# linear first order vectorial process: using a diagonal as the evolution matrix
function SP_1T_V(R_time::Array{Cdouble,1},L_source::Array{Cdouble,2},Nt::Int64, x0::Array{Cdouble,1}) # L_source could be made LowerTriangular or Diagonal or Cdouble
    function f_time(X::Array{Cdouble,1}) R_time.*X end
    SP_1T_V(f_time,L_source,Nt, x0)
end


# linear first order vectorial process: using an evolution model matrix
function SP_1T_V(R_time::Array{Cdouble,2},L_source::Array{Cdouble,2},Nt::Int64, x0::Array{Cdouble,1}) # L_source could be made LowerTriangular or Diagonal or Cdouble
    function f_time(X::Array{Cdouble,1}) R_time*X end
    SP_1T_V(f_time,L_source,Nt, x0)
end




# linear first order vectorial process asymptotically corrected: using a scalar as the evolution model (same evolution model for every elements)
function SP_1TC_V(r_time::Cdouble,L_proc::Array{Cdouble,2},Nt::Int64, x0::Array{Cdouble,1}) # L_source could be made LowerTriangular or Diagonal or Cdouble
    L_source = sqrt(1.0-r_time^2)*L_proc
    function f_time(X::Array{Cdouble,1}) r_time*X end
    SP_1T_V(f_time,L_source,Nt, x0)
end



# linear first order vectorial process asymptotically corrected: using a diagonal as the evolution matrix
function SP_1TC_V(R_time::Array{Cdouble,1},L_proc::Array{Cdouble,2},Nt::Int64, x0::Array{Cdouble,1}) # L_source could be made LowerTriangular or Diagonal or Cdouble
    L_source = chol_PS_L(L_proc*L_proc' - diagm(R_time)*L_proc*L_proc'*diagm(R_time))
    function f_time(X::Array{Cdouble,1}) R_time.*X end
    SP_1T_V(f_time,L_source,Nt, x0)
end





# linear first order vectorial process asymptotically corrected: using a diagonal as the evolution matrix
function SP_1TC_V(R_time::Array{Cdouble,2},L_proc::Array{Cdouble,2},Nt::Int64, x0::Array{Cdouble,1}) # L_source could be made LowerTriangular or Diagonal or Cdouble
    L_source = chol_PS_L(L_proc*L_proc' - R_time*L_proc*L_proc'*R_time')
    function f_time(X::Array{Cdouble,1}) R_time*X end
    SP_1T_V(f_time,L_source,Nt, x0)
end
