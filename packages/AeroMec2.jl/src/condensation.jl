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

# this file encodes the condensation mechanism
# it is assumed that the size dependence of the condensational growth rate is known and stored in ws.scale_GR
# the growth rate does not depend on the size (ws.GR could be a scalar instead of an array)
# The mechanism are described dimensionless: to retrieve the values with physical unit, one must multiply the results of CondensationGrowth! or jacobian_condensation! by ws.GR0

# condensation
function CondensationGrowth!(dx_cond::Array{Cdouble,1},ws::AeroSys,x::Array{Cdouble,1},zeta::Cdouble)
    ws.GR = (ws.GR0*ws.t0/ws.d0)*zeta # dimensionless growth rate
    # for a growth rate indenpend of the size
    dx_cond[1] = -ws.GR*ws.scale_GR[1]*x[1]
    # if ws.logS
        # dx_cond[2:end] = ws.GR*ws.scale_GR[2:end].*(ws.cst_r*x[1:end-1]-x[2:end])  # it a different type of discretization (growth estimated at the boundary of the bin range, it's not bad, but difficult to expalin)
        dx_cond[2:end] = ws.GR*(ws.scale_GR[1:end-1].*ws.cst_r.*x[1:end-1] - ws.scale_GR[2:end].*x[2:end])
    # else # cst_r = 1.0 if logS==false
    #     dx_cond[2:end] = ws.GR*(ws.scale_GR[1:end-1].*x[1:end-1] - ws.scale_GR[2:end].*x[2:end])
    # end
end

function CondensationGrowthAndEvap!(dx_cond::Array{Cdouble,1},ws::AeroSys,x::Array{Cdouble,1},zeta::Cdouble)
    ws.GR = (ws.GR0*ws.t0/ws.d0)*zeta # dimensionless growth rate
    # for a growth rate indenpend of the size
    if (zeta>=0)
        dx_cond[1] = -ws.GR*ws.scale_GR[1]*x[1]
        dx_cond[2:end] = ws.GR*(ws.scale_GR[1:end-1].*ws.cst_r.*x[1:end-1] - ws.scale_GR[2:end].*x[2:end])
    else
        dx_cond[end] = ws.GR*ws.scale_GR[end]*x[end]
        dx_cond[1:end-1] = ws.GR*(ws.scale_GR[1:end-1].*ws.cst_r.*x[1:end-1] - ws.scale_GR[2:end].*x[2:end])
    end
end

function CondensationGrowth!(dx_cond::Array{Cdouble,1},ws::AeroSys,x::Array{Cdouble,1},zeta::Array{Cdouble,1})
    GR_array = (ws.GR0*ws.t0/ws.d0)*zeta[:] # dimensionless growth rate
    ws.GR = GR_array[1]
    # for a growth rate indenpend of the size
    dx_cond[1] = -GR_array[1]*ws.scale_GR[1]*x[1]
    # dx_cond[2:end] = GR_array[2:end].*ws.scale_GR[2:end].*(ws.cst_r*x[1:end-1]-x[2:end]) # it a different type of discretization (growth estimated at the boundary of the bin range, it's not bad, but difficult to expalin)
    # if ws.logS
        dx_cond[2:end] = GR_array[1:end-1].*ws.scale_GR[1:end-1].*ws.cst_r.*x[1:end-1] - GR_array[2:end].*ws.scale_GR[2:end].*x[2:end]
    # else # cst_r = 1.0 if logS==false
    #     dx_cond[2:end] = GR_array[1:end-1].*ws.scale_GR[1:end-1].*x[1:end-1] - GR_array[2:end].*ws.scale_GR[2:end].*x[2:end]
    # end
end

function CondensationGrowthAndEvap!(dx_cond::Array{Cdouble,1},ws::AeroSys,x::Array{Cdouble,1},zeta::Array{Cdouble,1})
    GR_array = (ws.GR0*ws.t0/ws.d0)*zeta[:] # dimensionless growth rate
    ws.GR = GR_array[1]
    # positive growth rate (upwind)
    idx_pos = (GR_array.>=0.0)
    GR_pos = GR_array.*idx_pos;
    dx_cond[1] = -GR_pos[1]*ws.scale_GR[1]*x[1]
    dx_cond[2:end] = GR_pos[1:end-1].*ws.scale_GR[1:end-1].*ws.cst_r.*x[1:end-1] - GR_pos[2:end].*ws.scale_GR[2:end].*x[2:end]
    # negative growth rate (backwards)
    idx_neg = (GR_array.<0.0)
    GR_neg = GR_array.*idx_neg;
    dx_cond[end] = dx_cond[end] + GR_neg[end]*ws.scale_GR[end]*x[end]
    dx_cond[1:end-1] = dx_cond[1:end-1] + GR_neg[1:end-1].*ws.scale_GR[1:end-1].*ws.cst_r.*x[1:end-1] - GR_neg[2:end].*ws.scale_GR[2:end].*x[2:end]
end


# jacobian of the condensation term: OK
function jacobian_condensation!(F_ev_::Array{Cdouble,2},ws::AeroSys,x::Array{Cdouble,1},zeta::Cdouble) # assume that every element of the matrix is set to zero #LATER it might be improved by not computing the addition in the function
    ws.GR = (ws.GR0*ws.t0/ws.d0)*zeta # CGR(zeta) # dimensionless growth rate
    # for a growth rate indenpend of the size (known and stored in ws)
    F_ev_[1,1]   =     F_ev_[1,1]   - ws.GR*ws.scale_GR[1]
    for i in 2:ws.nbin
        # derivation w.r.t. the concentration distribution
        F_ev_[i,i]   = F_ev_[i,i]   - ws.GR*ws.scale_GR[i]
        # F_ev_[i,i-1] = F_ev_[i,i-1] + ws.GR*ws.scale_GR[i]*ws.cst_r
        # if ws.logS
            F_ev_[i,i-1] = F_ev_[i,i-1] + ws.GR*ws.scale_GR[i-1]*ws.cst_r
        # else # cst_r = 1.0 if logS==false
        #     F_ev_[i,i-1] = F_ev_[i,i-1] + ws.GR*ws.scale_GR[i-1]
        # end
    end
end

function jacobian_condensation!(F_ev_::Array{Cdouble,2},ws::AeroSys,x::Array{Cdouble,1},zeta::Array{Cdouble,1}) # assume that every element of the matrix is set to zero #LATER it might be improved by not computing the addition in the function
    GR_array = (ws.GR0*ws.t0/ws.d0)*zeta[:] # CGR(zeta) # dimensionless growth rate
    ws.GR = GR_array[1]
    # for a growth rate that denpends on the size
    F_ev_[1,1]   =     F_ev_[1,1]   - GR_array[1]*ws.scale_GR[1]
    for i in 2:ws.nbin
        # derivation w.r.t. the concentration distribution
        F_ev_[i,i]   = F_ev_[i,i]   - GR_array[i]*ws.scale_GR[i]
        # F_ev_[i,i-1] = F_ev_[i,i-1] + GR_array[i]*ws.scale_GR[i]*ws.cst_r # cf CondensationGrowth! comments
        # if ws.logS
            F_ev_[i,i-1] = F_ev_[i,i-1] + GR_array[i-1]*ws.scale_GR[i-1]*ws.cst_r
        # else # cst_r = 1.0 if logS==false
        #     F_ev_[i,i-1] = F_ev_[i,i-1] + GR_array[i-1]*ws.scale_GR[i-1]
        # end
    end
end


# error due to error in the parameters
"""
    cond_err(dt::Cdouble, e_g::Array{Cdouble,1},e_j::Cdouble,x::Array{Cdouble,1},delta::Array{Cdouble,1})

    computes the standard deviation of the error in the GDE, in terms of PSD, due to errors in the
    values of the growth rates and the nucleation rate.

    e_g   is an array of errors in the growth rates (e.g. sqrt(var(g-g^*)))
    e_j   is the value of the nucleation rate error (e.g. sqrt(var(J-J^*)))
    x     is the particle size distribution (number concentration in each bin)
    delta the widths of the bins of the size discretization
    dt    the time step

    Note: mind the units!
"""
function cond_err(dt::Cdouble, e_g::Array{Cdouble,1},e_j::Cdouble,x::Array{Cdouble,1},delta::Array{Cdouble,1})
    # e_g = (1.0e9*GR0*CGR(sqrt.(diag(myWSKF.o_fil[R_cond,R_cond]))))
    # e_j = Nucleation_rate(sqrt(myWSKF.o_fil[R_nuc_init,R_nuc_init]))
    # x = x0*myWSKF.x_fil[R_psd]
    # delta = 1.0e9delta
    # u_density.*e_g
    # (1.0e-9x0*myWSKF.x_fil[R_psd]./delta).*(1.0e9*GR0*CGR(sqrt.(diag(myWSKF.o_fil[R_cond,R_cond]))))
    # compute the size density from the size concentrations
    u_density = (x./delta)
    # approximation of the particle flux (the numerical flux cm^{-3} s^{-1}) at each bin boundary
    part_flux = [e_j ; u_density.*e_g]

    # return an approximation of sqrt(Var(dt*(part_flux[2:end] - part_flux[1:end-1])) assuming the elements of part_flux are not correlated, i.e. e_j and the entries of e_g are independents
    dt*sqrt.(part_flux[2:end].^2 + part_flux[1:end-1].^2)
end
