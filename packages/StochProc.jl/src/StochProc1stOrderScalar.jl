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


# non-linear first order process
function SP_1T(f_time::Function,s_source::Cdouble,Nt::Int64, x0::Cdouble=0.0)
    # allocate an array
    X = Array{Cdouble,1}(undef,Nt)
    # initi the time series
    X[1] = x0
    # generate the time series
    for i in 2:Nt
        X[i] = f_time(X[i-1]) + s_source*randn()
    end
    # return
    X
end

# linear first order process
function SP_1T(r_time::Cdouble,s_source::Cdouble,Nt::Int64, x0::Cdouble=0.0)
    function f_time(x::Cdouble) r_time*x end
    # return
    SP_1T(f_time,s_source,Nt, x0)
end

# linear first order process with corrected variance (not in time)
function SP_1TC(r_time::Cdouble,s_proc::Cdouble,Nt::Int64, x0::Cdouble=0.0)
    # linear function
    function f_time(x::Cdouble) r_time*x end
    # compute the amplitude of the source term
    s_source = s_proc*sqrt(1-r_time^2)
    # return
    SP_1T(f_time,s_source,Nt, x0)
end


# linear first order process with corrected variance and time correction
function SP_1TCV(r_time::Cdouble,s_proc::Cdouble,Nt::Int64, x0::Cdouble=0.0)
    # compute the amplitude of the source term
    s_source = s_proc*sqrt(1-r_time^2)
    # allocate an array
    X = Array{Cdouble,1}(undef,Nt)
    X_corr = Array{Cdouble,1}(undef,Nt)
    s_corr = Array{Cdouble,1}(undef,Nt)
    # initi the time series
    X[1] = x0
    s_time = s_source
    X_corr[1] = (s_proc/s_time)*x0
    s_corr[1] = s_proc/s_time
    # generate the time series
    for i in 2:Nt
        X[i] = r_time*X[i-1] + s_source*randn()
        s_time = sqrt(r_time^2*s_time^2 + s_source^2)
        s_corr[i] = s_proc/s_time
        X_corr[i] = s_corr[i]*X[i]
    end
    # return
    X_corr
end
