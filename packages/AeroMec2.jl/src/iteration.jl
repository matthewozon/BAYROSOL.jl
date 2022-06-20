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

# the equation being implemented:
# with dimension
# s = 1
# \frac{x_s}{\partial t} = -x_s \overset{N}{\underset{i=s+1}{\sum}} \beta_{s,i}x_i  - GR(t) \frac{S_s}{d_0} x_s - \gamma_s x_s + J
# and
# \forall s\in [\![2,N]\!]
# \frac{x_s}{\partial t} = -x_s \overset{N}{\underset{i=s+1}{\sum}} \beta_{s,i}x_i  + GR(t) \frac{S_s}{d_0} (r x_{s-1} - x_s) - \gamma_s x_s
# where x_s is the number concentration for the size bin s, N is the number of bins, \beta_{s,i} is the coagulation coefficient for the pair of particles size s and i, GR(t) is the growth rate due to condensing vapour at time t, S_s is the scaling factor of the condensation rate for bin s, d_0 a characteristic particle diameter, r is the constant ratio of the diameter series and \gamma_s the loss rate due to wall deposition.

"""
    **iter!**: the **!** in the syntax mean that the function may modify the arguments it is passed. It computes one time iteration of the aerosol system knowing the current state of the system (by definition, the state of the simulated system is the collection of concentrations, nucleation rate, coagulation coefficients, growth rates and loss rates). The time iteration is solved using a Euler discretization scheme.
"""
function iter!(dx_coa::Array{Cdouble,1},dx_con::Array{Cdouble,1},dx_nuc::Array{Cdouble,1},dx_los::Array{Cdouble,1},ws::AeroSys,x::Array{Cdouble,1},dt::Cdouble,zeta::Union{Cdouble,Array{Cdouble,1}},j::Cdouble,xi::Array{Cdouble,1})
    # coagulation
    if ws.is_coa
        if ws.is_coa_gain
            coagulation!(dx_coa,ws,x)
        else
            coagulation_loss!(dx_coa,ws,x)
        end
    else
        fill!(dx_coa,0.0)
    end
    # condensation
    if ws.is_con
        CondensationGrowth!(dx_con,ws,x,zeta)
    else
        fill!(dx_con,0.0)
    end
    # nucleation
    if ws.is_nuc
        Nucleation!(dx_nuc,ws,j)
    else
        fill!(dx_nuc,0.0)
    end
    # linear losses
    if ws.is_los
        WallDeposition!(dx_los,ws,x,xi)
    else
        fill!(dx_los,0.0)
    end

    # return the next size distribution
    x + (dt/ws.t0)*(dx_coa+dx_con+dx_nuc+dx_los)
end

"""
    **jacobian_GDE**: this function computes the Jacobian matrix **&delta;** of the discrete time evolution equation of the concentration N_i^{k+1} = f(N_i^k;J^k,\\beta^k,GR^k,\\gamma^k), i.e. the matrix whose entries are defined by: M_{i,j} = \\frac{\\partial N_i^{k+1}}{\\partial N_j^k}
"""
function jacobian_GDE(F_coa::Array{Cdouble,2},F_con::Array{Cdouble,2},F_los::Array{Cdouble,2},ws::AeroSys,x::Array{Cdouble,1},dt::Cdouble,zeta::Union{Cdouble,Array{Cdouble,1}},xi::Array{Cdouble,1})
    # coagulation
    if ws.is_coa
        if ws.is_coa_gain
            jacobian_coagulation!(F_coa,ws,x)
        else
            jacobian_coagulation_loss!(F_coa,ws,x)
        end
    else
        fill!(F_coa,0.0)
    end
    # condensation
    if ws.is_con
        jacobian_condensation!(F_con,ws,x,zeta)
    else
        fill!(F_con,0.0)
    end
    # linear losses
    if ws.is_los
        jacobian_wall_deposition!(F_los,ws,x,xi)
    else
        fill!(F_los,0.0)
    end

    # return de jacobian
    Matrix{Cdouble}(I,ws.nbin,ws.nbin) + (dt/ws.t0)*(F_coa+F_con+F_los)
end

"""
    discretization_err(dt::Cdouble,n_deriv::Array{Cdouble,1},g_star::Array{Cdouble,1},delta::Array{Cdouble,1})

    computes the standard deviation of the error in the GDE, in terms of PSD, due to the
    discretization in size and time (assuming all parameters known)

    n_deriv first order size-derivative of the size density
    g_star  growth rates
    delta   the widths of the bins of the size discretization
    dt      the time step

    Note: mind the units!
"""
function discretization_err(dt::Cdouble,n_deriv::Array{Cdouble,1},g_star::Array{Cdouble,1},delta::Array{Cdouble,1})
    # g_star = (GR0*CGR(myWSKF.x_fil[R_cond]))
    # delta = 1.0e9*delta
    cst_r = mean(delta[2:end]./delta[1:end-1])
    discretization_error = 0.5*(((sqrt(cst_r)-1)^2)/(cst_r-1.0))*(n_deriv.*g_star).*delta;
    dt*[abs(discretization_error[2]) + (1.0/cst_r)*abs(discretization_error[1]); abs.(discretization_error[2:end]) + (1.0/cst_r)*abs.(discretization_error[1:end-1])];
end
