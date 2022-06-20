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


# non-linear second order process
function SP_2T(f_time::Function,s_source::Cdouble,Nt::Int64, x0::Cdouble=0.0, x1::Cdouble=0.0)
    # allocate an array
    X = Array{Cdouble,1}(undef,Nt)
    # initi the time series
    X[1] = x0
    X[2] = x1
    # generate the time series
    for i in 3:Nt
        X[i] = f_time(X[i-1],X[i-2]) + s_source*randn()
    end
    # return
    X
end


# linear second order process with two distinct roots
function SP_2T(r_1::Cdouble,r_2::Cdouble,s_source::Cdouble,Nt::Int64, x0::Cdouble=0.0, x1::Cdouble=0.0)
    function f_time(x::Cdouble,y::Cdouble) x*(r_1+r_2) - y*r_1*r_2 end
    # return
    SP_2T(f_time,s_source,Nt, x0, x1)
end

# linear second order process with a double root
function SP_2T(r_time::Cdouble,s_source::Cdouble,Nt::Int64, x0::Cdouble=0.0, x1::Cdouble=0.0)
    SP_2T(r_time,r_time,s_source,Nt, x0, x1)
end


# linear second order process with asymptotic variance correction
function SP_2TC(r_1::Cdouble,r_2::Cdouble,s_proc::Cdouble,Nt::Int64, x0::Cdouble=0.0, x1::Cdouble=0.0)
    # compute the amplitude of the source term
    s_source = s_proc*sqrt(1.0 -(r_1+r_2)^2 -(r_1*r_2)^2 +2.0*r_1*r_2*((r_1+r_2)^2)/(1.0+r_1*r_2) )
    # return
    SP_2T(r_1,r_2,s_source,Nt, x0, x1)
end

function SP_2TC(r_time::Cdouble,s_proc::Cdouble,Nt::Int64, x0::Cdouble=0.0, x1::Cdouble=0.0)
    # return
    SP_2TC(r_time,r_time,s_proc,Nt, x0, x1)
end

# time correction using the variance: it conserves the smoothness of the process
function SP_2TCV(r_1::Cdouble,r_2::Cdouble,s_proc::Cdouble,Nt::Int64, x0::Cdouble=0.0, x1::Cdouble=0.0)
    # time evolution model
    B = [r_1+r_2  -r_1*r_2; 1.0 0.0]
    # compute the amplitude of the source term
    s_source = s_proc*sqrt(1.0 -(r_1+r_2)^2 -(r_1*r_2)^2 +2.0*r_1*r_2*((r_1+r_2)^2)/(1.0+r_1*r_2) )
    G_source = [s_source^2 0.0; 0.0 0.0]
    L_source = [s_source 0.0; 0.0 0.0]
    # allocate some arrays
    X      = Array{Cdouble,2}(undef,2,Nt)
    X_corr = Array{Cdouble,2}(undef,2,Nt)
    G_cov  = Array{Cdouble,3}(undef,2,2,Nt)
    # init the time series
    X[1,1] = x1
    X[2,1] = x0
    G_cov[1,1,1] = s_proc^2; G_cov[1,2,1] = 0.0
    G_cov[2,1,1] = 0.0;      G_cov[2,2,1] = 0.0
    X_corr[1,1] = x1
    X_corr[2,1] = x0
    # generate the time series
    for i in 2:Nt
        # propagate the state
        X[:,i] = B*X[:,i-1] + L_source*randn(2)
        # propagate the covariance
        G_cov[:,:,i] = B*G_cov[:,:,i-1]*B' + G_source
        # correction
        X_corr[:,i] = s_proc*X[:,i]./sqrt.([G_cov[1,1,i]; G_cov[2,2,i]])
    end
    # return
    X_corr[1,:]
end


# time correction using the covariance matrices of the first order vectorial system (for the modelization of the second order): it changes the smoothness of the process
function SP_2TCV_mat(r_1::Cdouble,r_2::Cdouble,s_proc::Cdouble,Nt::Int64, x0::Cdouble=0.0, x1::Cdouble=0.0)
    # time evolution model
    B = [r_1+r_2  -r_1*r_2; 1.0 0.0]
    # compute the amplitude of the source term
    s_source = s_proc*sqrt(1.0 -(r_1+r_2)^2 -(r_1*r_2)^2 +2.0*r_1*r_2*((r_1+r_2)^2)/(1.0+r_1*r_2) )
    G_source = [s_source^2 0.0; 0.0 0.0]
    L_source = [s_source 0.0; 0.0 0.0]
    # desired covariance of the process
    G_infty = (s_proc^2)*[1.0 r_1*r_2/(1.0 + r_1*r_2); r_1*r_2/(1.0 + r_1*r_2) 1.0]
    # A_infty = sqrtm(G_infty)
    A_infty = sqrt(G_infty)
    # allocate some arrays
    X      = Array{Cdouble,2}(undef,2,Nt)
    X_corr = Array{Cdouble,2}(undef,2,Nt)
    G_cov  = Array{Cdouble,3}(undef,2,2,Nt)
    # init the time series
    X[1,1] = x1
    X[2,1] = x0
    G_cov[1,1,1] = s_proc^2; G_cov[1,2,1] = 0.0
    G_cov[2,1,1] = 0.0;      G_cov[2,2,1] = 0.0
    X_corr[1,1] = x1
    X_corr[2,1] = x0
    # generate the time series
    for i in 2:Nt
        # propagate the state
        X[:,i] = B*X[:,i-1] + L_source*randn(2)
        # propagate the covariance
        G_cov[:,:,i] = B*G_cov[:,:,i-1]*B' + G_source
        # correction
        # X_corr[:,i] = A_infty*inv(sqrtm(G_cov[:,:,i]))*X[:,i]
        X_corr[:,i] = A_infty*inv(sqrt(G_cov[:,:,i]))*X[:,i]
    end
    # return
    X_corr[1,:]
end


function SP_2TCV(r_time::Cdouble,s_proc::Cdouble,Nt::Int64, x0::Cdouble=0.0, x1::Cdouble=0.0)
    SP_2TCV(r_time,r_time,s_proc,Nt, x0, x1)
end
